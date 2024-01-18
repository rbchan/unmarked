get_vars <- function(inp){
  if(is.list(inp)){
    out <- unique(unlist(lapply(inp, all.vars)))
  } else {
    out <- all.vars(inp)
  }
  names(out) <- out
  out
}

var_data <- function(var, guide, n){
  out <- rep(NA, n)
  gv <- guide[[var]]
  if(is.null(gv)){
    out <- stats::rnorm(n, 0, 1)
  } else if(inherits(gv, "factor")){
    levs <- levels(gv)
    out <- factor(sample(levs, n, replace=TRUE), levels=levs)
  } else{
    gv$n <- n
    out <- do.call(gv$dist, gv[!names(gv)=="dist"])
  }
  out
}

generate_data <- function(formulas, guide, n){
  vars <- get_vars(formulas)
  if(length(vars)==0) return(NULL)
  as.data.frame(lapply(vars, var_data, guide=guide, n=n))
}

capitalize <- function(inp){
  paste0(toupper(substring(inp,1,1)),
         substring(inp,2,nchar(inp)))
}

parse_func_name <- function(inp){
  if(!is.character(inp)){
    stop("Argument must be a character string", call.=FALSE)
  }
  capitalize(inp)
}

blank_umFit <- function(fit_function){
  type <- parse_func_name(fit_function)
  type <- ifelse(type=="Pcount", "PCount", type)
  type <- ifelse(type=="MultinomPois", "MPois", type)
  type <- ifelse(type=="Distsamp", "DS", type)
  type <- ifelse(type=="Colext", "ColExt", type)
  type <- ifelse(type=="Gdistsamp", "GDS", type)
  type <- ifelse(type=="Gpcount", "GPC", type)
  type <- ifelse(type=="Gmultmix", "GMM", type)
  type <- ifelse(type=="PcountOpen", "PCO", type)
  type <- ifelse(type=="DistsampOpen", "DSO", type)
  type <- ifelse(type=="MultmixOpen", "MMO", type)
  type <- ifelse(type=="Gdistremoval", "GDR", type)
  type <- paste0("unmarkedFit", type)
  new(type)
}


setMethod("simulate", "character",
  function(object, nsim=1, seed=NULL, formulas, coefs=NULL, design, guide=NULL, ...){
  model <- blank_umFit(object)
  fit <- suppressWarnings(simulate_fit(model, formulas, guide, design, ...))
  coefs <- check_coefs(coefs, fit)
  #fit <- replace_sigma(coefs, fit)
  coefs <- generate_random_effects(coefs, fit)
  fit <- replace_estimates(fit, coefs)
  ysims <- suppressWarnings(simulate(fit, nsim))
  umf <- fit@data
  # fix this
  umfs <- lapply(ysims, function(x){
    if(object=="occuMulti"){
      umf@ylist <- x
    } else if(object=="gdistremoval"){
      umf@yDistance=x$yDistance
      umf@yRemoval=x$yRemoval
    } else if(object == "IDS"){
       out <- list()
       out$ds <- fit@data
       out$ds@y <- x$ds
       if("pc" %in% names(fit)){
         out$pc <- fit@dataPC
         out$pc@y <- x$pc
       }
       if("oc" %in% names(fit)){
         out$oc <- fit@dataOC
         out$oc@y <- x$oc
       }
       umf <- out
    } else {
      umf@y <- x
    }
    umf
  })
  if(length(umfs)==1) umfs <- umfs[[1]]
  umfs
})

# Insert specified random effects SD into proper S4 slot in model object
# This is mostly needed by GDR which uses the SD to calculate
# N with E_loglam (this is currently disabled so the function is not needed)
#replace_sigma <- function(coefs, fit){
#  required_subs <- names(fit@estimates@estimates)
#  formulas <- sapply(names(fit), function(x) get_formula(fit, x))
#  rand <- lapply(formulas, lme4::findbars)
#  if(!all(sapply(rand, is.null))){
#    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
#    for (i in required_subs){
#      if(!is.null(rand[[i]][[1]])){
#        signame <- rvar[[i]]
#        old_coefs <- coefs[[i]]
#        fit@estimates@estimates[[i]]@randomVarInfo$estimates <- coefs[[i]][[signame]]
#      }
#    }
#  }
#  fit
#}

generate_random_effects <- function(coefs, fit){
  required_subs <- names(fit@estimates@estimates)
  formulas <- sapply(names(fit), function(x) get_formula(fit, x))
  rand <- lapply(formulas, lme4::findbars)
  if(!all(sapply(rand, is.null))){
    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
    for (i in required_subs){
      if(!is.null(rand[[i]][[1]])){
        signame <- rvar[[i]]
        old_coefs <- coefs[[i]]
        new_coefs <- old_coefs[names(old_coefs)!=signame]

        # Find levels of factor variable
        if(signame %in% names(siteCovs(fit@data))){
          lvldata <- siteCovs(fit@data)[[signame]]
        } else if(signame %in% names(obsCovs(fit@data))){
          lvldata <- obsCovs(fit@data)[[signame]]
        } else if(methods::.hasSlot(fit@data, "yearlySiteCovs") && signame %in% names(yearlySiteCovs(fit@data))){
          lvldata <- yearlySiteCovs(fit@data)[[signame]]
        } else {
          stop("Random effect covariate missing from data", call.=FALSE)
        }

        if(!is.factor(lvldata)){
          stop("Random effect covariates must be specified as factors with guide argument", call.=FALSE)
        }
        b <- stats::rnorm(length(levels(lvldata)), 0, old_coefs[signame])
        names(b) <- rep(paste0("b_",i), length(b))
        new_coefs <- c(new_coefs, b)
        coefs[[i]] <- new_coefs
      }
    }
  }
  coefs
}


setGeneric("get_umf_components", function(object, ...) standardGeneric("get_umf_components"))

setMethod("get_umf_components", "unmarkedFit",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$state, guide, design$M)
  oc <- generate_data(formulas$det, guide, design$J*design$M)
  yblank <- matrix(0, design$M, design$J)
  list(y=yblank, siteCovs=sc, obsCovs=oc)
})


setGeneric("simulate_fit", function(object, ...) standardGeneric("simulate_fit"))

setMethod("simulate_fit", "unmarkedFitOccu",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFrameOccu(y=parts$y, siteCovs=parts$siteCovs,
                           obsCovs=parts$obsCovs)
  occu(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
       data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("simulate_fit", "unmarkedFitPCount",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFramePCount(y=parts$y, siteCovs=parts$siteCovs,
                             obsCovs=parts$obsCovs)
  args <- list(...)
  K <- ifelse(is.null(args$K), 100, args$K)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)
  pcount(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
       data=umf, mixture=mixture, K=K, se=FALSE, control=list(maxit=1))
})

setMethod("simulate_fit", "unmarkedFitOccuRN",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedFrameOccu(y=parts$y, siteCovs=parts$siteCovs,
                           obsCovs=parts$obsCovs)
  occuRN(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
       data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("get_umf_components", "unmarkedFitMPois",
          function(object, formulas, guide, design, ...){
  args <- list(...)
  sc <- generate_data(formulas$state, guide, design$M)
  oc <- generate_data(formulas$det, guide, design$J*design$M)
  J <- ifelse(args$type=="double", 3, design$J)
  yblank <- matrix(0, design$M, design$J)
  list(y=yblank, siteCovs=sc, obsCovs=oc)
})

setMethod("simulate_fit", "unmarkedFitMPois",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  type <- ifelse(is.null(args$type), "removal", args$type)
  umf <- unmarkedFrameMPois(y=parts$y, siteCovs=parts$siteCovs,
                           obsCovs=parts$obsCovs, type=type)
  multinomPois(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
               data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("get_umf_components", "unmarkedFitDS",
          function(object, formulas, guide, design, ...){
  #args <- list(...)
  sc <- generate_data(formulas$state, guide, design$M)
  sc2 <- generate_data(formulas$det, guide, design$M)
  dat <- list(sc, sc2)
  keep <- sapply(dat, function(x) !is.null(x))
  dat <- dat[keep]
  sc <- do.call(cbind, dat)
  yblank <- matrix(0, design$M, design$J)
  list(y=yblank, siteCovs=sc)
})

setMethod("simulate_fit", "unmarkedFitDS",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(is.null(args$tlength)) args$tlength <- 0
  umf <- unmarkedFrameDS(y=parts$y, siteCovs=parts$siteCovs,
                         tlength=args$tlength, survey=args$survey, unitsIn=args$unitsIn,
                         dist.breaks=args$dist.breaks)
  keyfun <- ifelse(is.null(args$keyfun), "halfnorm", args$keyfun)
  output <- ifelse(is.null(args$output), "density", args$output)
  unitsOut <- ifelse(is.null(args$unitsOut), "ha", args$unitsOut)

  distsamp(as.formula(paste(deparse(formulas$det), deparse(formulas$state))),
               data=umf, se=FALSE, control=list(maxit=1), keyfun=keyfun,
          output=output, unitsOut=unitsOut)
})


setMethod("get_umf_components", "unmarkedFitColExt",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$psi, guide, design$M)
  ysc <- generate_data(list(formulas$col, formulas$ext), guide, design$M*design$T)
  oc <- generate_data(formulas$det, guide, design$J*design$M*design$T)
  yblank <- matrix(0, design$M, design$T*design$J)
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc, obsCovs=oc)
})


setMethod("simulate_fit", "unmarkedFitColExt",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  umf <- unmarkedMultFrame(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           obsCovs=parts$obsCovs, numPrimary=design$T)
  colext(psiformula=formulas$psi, gammaformula=formulas$col,
         epsilonformula=formulas$ext,pformula=formulas$det,
         data=umf, se=FALSE, control=list(maxit=1))
})

setMethod("get_umf_components", "unmarkedFitOccuTTD",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$psi, guide, design$M)
  ysc <- NULL
  if(design$T>1){
    ysc <- generate_data(list(formulas$col, formulas$ext), guide, design$M*design$T)
  }
  oc <- generate_data(formulas$det, guide, design$J*design$M*design$T)
  yblank <- matrix(0, design$M, design$T*design$J)
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc, obsCovs=oc)
})


setMethod("simulate_fit", "unmarkedFitOccuTTD",
  function(object, formulas, guide, design, ...){
  if(is.null(design$T)) design$T <- 1
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  umf <- unmarkedFrameOccuTTD(y=parts$y,
                              surveyLength=args$surveyLength,
                              siteCovs=parts$siteCovs,
                              yearlySiteCovs=parts$yearlySiteCovs,
                              obsCovs=parts$obsCovs, numPrimary=design$T)
  linkPsi <- ifelse(is.null(args$linkPsi), "logit", args$linkPsi)
  ttdDist <- ifelse(is.null(args$ttdDist), "exp", args$ttdDist)
  occuTTD(psiformula=formulas$psi, gammaformula=formulas$col,
         epsilonformula=formulas$ext,detformula=formulas$det,
         linkPsi=linkPsi, ttdDist=ttdDist,
         data=umf, se=FALSE, control=list(maxit=1))
})


setMethod("get_umf_components", "unmarkedFitGMM",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$lambda, guide, design$M)
  ysc <- generate_data(formulas$phi, guide, design$M*design$T)
  yblank <- matrix(0, design$M, design$T*design$J)
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc)
})


setMethod("simulate_fit", "unmarkedFitGDS",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(args$survey=="line"){
    umf <- unmarkedFrameGDS(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           numPrimary=design$T,
                           tlength=args$tlength, survey=args$survey,
                           unitsIn=args$unitsIn, dist.breaks=args$dist.breaks)
  } else if(args$survey=="point"){
    umf <- unmarkedFrameGDS(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           numPrimary=design$T, survey=args$survey,
                           unitsIn=args$unitsIn, dist.breaks=args$dist.breaks)
  }

  keyfun <- ifelse(is.null(args$keyfun), "halfnorm", args$keyfun)
  output <- ifelse(is.null(args$output), "density", args$output)
  unitsOut <- ifelse(is.null(args$unitsOut), "ha", args$unitsOut)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)
  K <- ifelse(is.null(args$K), 100, args$K)

  gdistsamp(lambdaformula=formulas$lambda, phiformula=formulas$phi,
            pformula=formulas$det, data=umf, keyfun=keyfun, output=output,
            unitsOut=unitsOut, mixture=mixture, K=K,
            se=FALSE, control=list(maxit=1))
})

setMethod("simulate_fit", "unmarkedFitGPC",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  umf <- unmarkedFrameGPC(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           numPrimary=design$T)
  K <- ifelse(is.null(args$K), 100, args$K)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)

  gpcount(lambdaformula=formulas$lambda, phiformula=formulas$phi,
          pformula=formulas$det, data=umf, mixture=mixture, K=K,
          se=FALSE, control=list(maxit=1))
})


setMethod("simulate_fit", "unmarkedFitGMM",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  umf <- unmarkedFrameGMM(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           numPrimary=design$T, type=args$type)
  K <- ifelse(is.null(args$K), 100, args$K)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)

  gmultmix(lambdaformula=formulas$lambda, phiformula=formulas$phi,
           pformula=formulas$det, data=umf, mixture=mixture, K=K,
           se=FALSE, control=list(maxit=1))
})


setMethod("get_umf_components", "unmarkedFitDailMadsen",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$lambda, guide, design$M)
  ysc <- generate_data(list(formulas$gamma, formulas$omega), guide, design$M*design$T)
  oc <- generate_data(formulas$det, guide, design$M*design$T*design$J)
  yblank <- matrix(0, design$M, design$T*design$J)
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc, obsCovs=oc)
})

setMethod("simulate_fit", "unmarkedFitPCO",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(is.null(args$primaryPeriod)){
    args$primaryPeriod <- matrix(1:design$T, design$M, design$T, byrow=TRUE)
  }
  umf <- unmarkedFramePCO(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           numPrimary=design$T, primaryPeriod=args$primaryPeriod)
  K <- ifelse(is.null(args$K), 100, args$K)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)
  dynamics <- ifelse(is.null(args$dynamics), "constant", args$dynamics)
  fix <- ifelse(is.null(args$fix), "none", args$fix)
  immigration <- ifelse(is.null(args$immigration), FALSE, args$immigration)
  iotaformula <- args$iotaformula
  if(is.null(iotaformula)) iotaformula <- ~1

  pcountOpen(lambdaformula=formulas$lambda, gammaformula=formulas$gamma,
             omegaformula=formulas$omega, pformula=formulas$det,
             data=umf, mixture=mixture, K=K, dynamics=dynamics, fix=fix,
             se=FALSE, method='SANN', control=list(maxit=1), immigration=immigration,
             iotaformula=iotaformula)
})

setMethod("simulate_fit", "unmarkedFitMMO",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(is.null(args$primaryPeriod)){
    args$primaryPeriod <- matrix(1:design$T, design$M, design$T, byrow=TRUE)
  }
  umf <- unmarkedFrameMMO(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           type=args$type,
                           numPrimary=design$T, primaryPeriod=args$primaryPeriod)
  K <- ifelse(is.null(args$K), 100, args$K)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)
  dynamics <- ifelse(is.null(args$dynamics), "constant", args$dynamics)
  fix <- ifelse(is.null(args$fix), "none", args$fix)
  immigration <- ifelse(is.null(args$immigration), FALSE, args$immigration)
  iotaformula <- args$iotaformula
  if(is.null(iotaformula)) iotaformula <- ~1

  multmixOpen(lambdaformula=formulas$lambda, gammaformula=formulas$gamma,
             omegaformula=formulas$omega, pformula=formulas$det,
             data=umf, mixture=mixture, K=K, dynamics=dynamics, fix=fix,
             se=FALSE, method='SANN', control=list(maxit=1), immigration=immigration,
             iotaformula=iotaformula)
})

setMethod("get_umf_components", "unmarkedFitDSO",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(formulas$lambda, guide, design$M)
  ysc <- generate_data(list(formulas$gamma, formulas$omega, formulas$det),
                       guide, design$M*design$T)
  yblank <- matrix(0, design$M, design$T*design$J)
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc)
})

setMethod("simulate_fit", "unmarkedFitDSO",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(is.null(args$primaryPeriod)){
    args$primaryPeriod <- matrix(1:design$T, design$M, design$T, byrow=TRUE)
  }
  if(args$survey=="line"){
  umf <- unmarkedFrameDSO(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,
                           tlength=args$tlength, survey=args$survey,
                           unitsIn=args$unitsIn, dist.breaks=args$dist.breaks,
                           numPrimary=design$T, primaryPeriod=args$primaryPeriod)
  } else if(args$survey == "point"){
  umf <- unmarkedFrameDSO(y=parts$y, siteCovs=parts$siteCovs,
                           yearlySiteCovs=parts$yearlySiteCovs,survey=args$survey,
                           unitsIn=args$unitsIn, dist.breaks=args$dist.breaks,
                           numPrimary=design$T, primaryPeriod=args$primaryPeriod)
  }
  K <- ifelse(is.null(args$K), 100, args$K)
  keyfun <- ifelse(is.null(args$keyfun), "halfnorm", args$keyfun)
  output <- ifelse(is.null(args$output), "density", args$output)
  unitsOut <- ifelse(is.null(args$unitsOut), "ha", args$unitsOut)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)
  dynamics <- ifelse(is.null(args$dynamics), "constant", args$dynamics)
  fix <- ifelse(is.null(args$fix), "none", args$fix)
  immigration <- ifelse(is.null(args$immigration), FALSE, args$immigration)
  iotaformula <- args$iotaformula
  if(is.null(iotaformula)) iotaformula <- ~1
  distsampOpen(lambdaformula=formulas$lambda, gammaformula=formulas$gamma,
             omegaformula=formulas$omega, pformula=formulas$det,
             keyfun=keyfun, unitsOut=unitsOut, output=output,
             data=umf, mixture=mixture, K=K, dynamics=dynamics, fix=fix,
             se=FALSE, method='SANN', control=list(maxit=1), immigration=immigration,
             iotaformula=iotaformula)
})


setMethod("get_umf_components", "unmarkedFitOccuMulti",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(lapply(formulas$state, as.formula), guide, design$M)
  oc <- generate_data(lapply(formulas$det, as.formula), guide, design$J*design$M)
  nspecies <- length(formulas$det)
  yblank <- lapply(1:nspecies, function(x) matrix(0, design$M, design$J))
  list(y=yblank, siteCovs=sc, obsCovs=oc)
})

setMethod("simulate_fit", "unmarkedFitOccuMulti",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  if(is.null(args$maxOrder)) args$maxOrder <- length(parts$y)
  umf <- unmarkedFrameOccuMulti(y=parts$y, siteCovs=parts$siteCovs,
                                obsCovs=parts$obsCovs, maxOrder=args$maxOrder)
  occuMulti(formulas$det, formulas$state, data=umf, maxOrder=args$maxOrder,
            se=FALSE, control=list(maxit=1))
})

setMethod("get_umf_components", "unmarkedFitOccuMS",
          function(object, formulas, guide, design, ...){
  sc <- generate_data(lapply(formulas$state, as.formula), guide, design$M)
  ysc <- NULL
  if(!is.null(formulas$phi)){
    ysc <- generate_data(lapply(formulas$phi, as.formula), guide, design$M*design$T*design$J)
  }
  oc <- generate_data(lapply(formulas$det, as.formula), guide, design$J*design$M)
  nspecies <- length(formulas$det)
  yblank <- matrix(0, design$M, design$T*design$J)
  yblank[1,1] <- 2 # To bypass sanity checker in unmarkedFrameOccuMS
  list(y=yblank, siteCovs=sc, yearlySiteCovs=ysc, obsCovs=oc)
})

setMethod("simulate_fit", "unmarkedFitOccuMS",
  function(object, formulas, guide, design, ...){
  if(is.null(design$T)) design$T <- 1
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  umf <- unmarkedFrameOccuMS(y=parts$y, siteCovs=parts$siteCovs,
                                yearlySiteCovs=parts$yearlySiteCovs,
                                obsCovs=parts$obsCovs, numPrimary=design$T)
  if(is.null(args$parameterization)) args$parameterization <- "multinomial"
  occuMS(formulas$det, formulas$state, formulas$phi, data=umf,
         parameterization=args$parameterization,
         se=FALSE, control=list(maxit=1))
})

setMethod("get_umf_components", "unmarkedFitGDR",
          function(object, formulas, guide, design, ...){
  if(any(! c("M","Jdist","Jrem") %in% names(design))){
    stop("Required design components are M, Jdist, and Jrem")
  }
  sc <- generate_data(list(formulas$lambda, formulas$dist), guide, design$M)
  ysc <- NULL
  if(design$T > 1){
    ysc <- generate_data(formulas$phi, guide, design$M*design$T)
  }
  oc <- generate_data(formulas$rem, guide, design$M*design$T*design$Jrem)

  list(yDistance=matrix(0, design$M, design$T*design$Jdist),
       yRemoval=matrix(0, design$M, design$T*design$Jrem),
       siteCovs=sc, yearlySiteCovs=ysc, obsCovs=oc)
})

setMethod("simulate_fit", "unmarkedFitGDR",
  function(object, formulas, guide, design, ...){
  if(is.null(design$T)) design$T <- 1
  if(is.null(formulas$phi)) formulas$phi <- ~1
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)
  umf <- unmarkedFrameGDR(yDistance=parts$yDistance, yRemoval=parts$yRemoval,
                          numPrimary=design$T, siteCovs=parts$siteCovs,
                          obsCovs=parts$obsCovs, yearlySiteCovs=parts$yearlySiteCovs,
                          dist.breaks=args$dist.breaks, unitsIn=args$unitsIn,
                          period.lengths=args$period.lengths)

  keyfun <- ifelse(is.null(args$keyfun), "halfnorm", args$keyfun)
  output <- ifelse(is.null(args$output), "density", args$output)
  unitsOut <- ifelse(is.null(args$unitsOut), "ha", args$unitsOut)
  mixture <- ifelse(is.null(args$mixture), "P", args$mixture)
  K <- ifelse(is.null(args$K), 100, args$K)

  gdistremoval(lambdaformula=formulas$lambda, phiformula=formulas$phi,
               removalformula=formulas$rem, distanceformula=formulas$dist,
               data=umf, keyfun=keyfun, output=output, unitsOut=unitsOut,
               mixture=mixture, K=K, se=FALSE, control=list(maxit=1), method='L-BFGS-B')
})

# For simulating entirely new datasets
setMethod("get_umf_components", "unmarkedFitIDS",
          function(object, formulas, guide, design, ...){

  # Distance sampling dataset
  sc_ds_lam <- generate_data(formulas$lam, guide, design$Mds)
  sc_ds_det <- generate_data(formulas$ds, guide, design$Mds)
  dat_ds <- list(sc_ds_lam, sc_ds_det)
  if(!is.null(formulas$phi)){
    sc_ds_phi <- generate_data(formulas$phi, guide, design$Mds)
    dat_ds <- c(dat_ds, list(sc_ds_phi))
  }
  keep <- sapply(dat_ds, function(x) !is.null(x))
  dat_ds <- dat_ds[keep]
  sc_ds <- do.call(cbind, dat_ds)
  yblank_ds <- matrix(1, design$Mds, design$J)

  # Point count dataset
  sc_pc <- yblank_pc <- NULL
  if(!is.null(design$Mpc) && design$Mpc > 0){
    if(is.null(formulas$pc)) form_pc <- formulas$ds
    sc_pc_lam <- generate_data(formulas$lam, guide, design$Mpc)
    sc_pc_det <- generate_data(form_pc, guide, design$Mpc)
    sc_pc <- list(sc_pc_lam, sc_pc_det)
    if(!is.null(formulas$phi)){
      sc_pc_phi <- generate_data(formulas$phi, guide, design$Mpc)
      sc_pc <- c(sc_pc, list(sc_pc_phi))
    }
    keep <- sapply(sc_pc, function(x) !is.null(x))
    sc_pc <- sc_pc[keep]
    sc_pc <- do.call(cbind, sc_pc)
    yblank_pc <- matrix(1, design$Mpc, 1)
  }

  # Presence/absence dataset
  sc_oc <- yblank_oc <- NULL
  if(!is.null(design$Moc) && design$Moc > 0){
    if(is.null(formulas$oc)){
      form_oc <- formulas$ds
    } else {
      form_oc <- formulas$oc
    }
    sc_oc_lam <- generate_data(formulas$lam, guide, design$Moc)
    sc_oc_det <- generate_data(form_oc, guide, design$Moc)
    sc_oc <- list(sc_oc_lam, sc_oc_det)
    if(!is.null(formulas$phi)){
      sc_oc_phi <- generate_data(formulas$phi, guide, design$Moc)
      sc_oc <- c(sc_oc, list(sc_oc_phi))
    }
    keep <- sapply(sc_oc, function(x) !is.null(x))
    sc_oc <- sc_oc[keep]
    sc_oc <- do.call(cbind, sc_oc)
    yblank_oc <- matrix(1, design$Moc, 1)
  }

  mget(c("yblank_ds", "sc_ds", "yblank_pc", "sc_pc", "yblank_oc", "sc_oc"))
})


setMethod("simulate_fit", "unmarkedFitIDS",
  function(object, formulas, guide, design, ...){
  parts <- get_umf_components(object, formulas, guide, design, ...)
  args <- list(...)

  args$tlength <- 0
  args$survey <- "point"

  # Distance sampling dataset
  umf_ds <- unmarkedFrameDS(y=parts$yblank_ds, siteCovs=parts$sc_ds,
                            tlength=args$tlength, survey=args$survey,
                            unitsIn=args$unitsIn,
                            dist.breaks=args$dist.breaks)
  # Point count dataset
  umf_pc <- NULL
  if(!is.null(design$Mpc) && design$Mpc > 0){
    umf_pc <- unmarkedFramePCount(y=parts$yblank_pc, siteCovs=parts$sc_pc)
  }

  # Occupancy dataset
  umf_oc <- NULL
  if(!is.null(design$Moc) && design$Moc > 0){
    umf_oc <- unmarkedFrameOccu(y=parts$yblank_oc, siteCovs=parts$sc_oc)
  }

  keyfun <- ifelse(is.null(args$keyfun), "halfnorm", args$keyfun)
  unitsOut <- ifelse(is.null(args$unitsOut), "ha", args$unitsOut)
  K <- ifelse(is.null(args$K), 300, args$K)
  if(is.null(args$maxDistPC)) args$maxDistPC <- max(args$dist.breaks)
  if(is.null(args$maxDistOC)) args$maxDistOC <- max(args$dist.breaks)

  IDS(lambdaformula = formulas$lam,
      detformulaDS = formulas$ds,
      detformulaPC = formulas$pc, detformulaOC = formulas$oc,
      dataDS = umf_ds, dataPC = umf_pc, dataOC = umf_oc,
      availformula = formulas$phi,
      durationDS = args$durationDS, durationPC = args$durationPC,
      durationOC = args$durationOC,
      maxDistPC = args$maxDistPC, maxDistOC = args$maxDistOC,
      keyfun=keyfun, unitsOut=unitsOut, K=K ,control=list(maxit=1))
})



