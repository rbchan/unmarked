setClass("unmarkedPower",
  representation(call="call", data="unmarkedFrame", M="numeric",
                 J="numeric", T="numeric", coefs="list", estimates="list",
                 alpha="numeric", nulls="list")
)

powerAnalysis <- function(object, coefs=NULL, design=NULL, alpha=0.05, nulls=list(),
                          datalist=NULL,
                          nsim=ifelse(is.null(datalist), 100, length(datalist)),
                          parallel=FALSE){

  stopifnot(inherits(object, "unmarkedFit"))

  submodels <- names(object@estimates@estimates)
  coefs <- check_coefs(coefs, object)
  coefs <- generate_random_effects(coefs, object)
  fit_temp <- replace_estimates(object, coefs)

  T <- 1
  bdata <- NULL
  if(!is.null(datalist)){
    if(length(datalist) != nsim){
      stop("Length of data list must equal value of nsim", call.=FALSE)
    }
    tryCatch({test <- update(object, data=datalist[[1]], se=FALSE,
                        control=list(maxit=1))
    }, error=function(e){
      stop("Incorrect format of entries in datalist", call.=FALSE)
    })
    bdata <- datalist
    M <- numSites(bdata[[1]])
    sims <- lapply(bdata, function(x){
      #fit_temp@data <- x
      #temporary workaround - not necessary??
      #if(methods::.hasSlot(fit_temp, "knownOcc")){
      #  fit_temp@knownOcc <- rep(FALSE, M)
      #}
      #simulate(fit_temp, 1)[[1]]
      if(inherits(x, "unmarkedFrameOccuMulti")){
        return(x@ylist)
      } else if(inherits(x, "unmarkedFrameGDR")){
        return(list(yDistance=x@yDistance, yRemoval=x@yRemoval))
      } else {
        return(x@y)
      }
    })
    if(methods::.hasSlot(bdata[[1]], "numPrimary")){
      T <- bdata[[1]]@numPrimary
    }
    J <- obsNum(bdata[[1]]) / T
  } else if(is.null(design)){
      sims <- simulate(fit_temp, nsim)
      M <- numSites(object@data)
      if(methods::.hasSlot(object@data, "numPrimary")){
        T <- object@data@numPrimary
      }
      J <- obsNum(object@data) / T
  } else {
    bdata <- bootstrap_data(fit_temp@data, nsim, design)
    sims <- lapply(bdata, function(x){
      fit_temp@data <- x
      #temporary workaround
      if(methods::.hasSlot(fit_temp, "knownOcc")){
        fit_temp@knownOcc <- rep(FALSE, design$M)
      }
      simulate(fit_temp, 1)[[1]]
    })
    M <- design$M
    if(methods::.hasSlot(fit_temp@data, "numPrimary")){
      T <- fit_temp@data@numPrimary
    }
    J <- design$J
  }

  cl <- NULL
  if(parallel){
    cl <- parallel::makeCluster(parallel::detectCores()-1)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(unmarked))
  }

  if(!is.null(options()$unmarked_shiny)&&options()$unmarked_shiny){
    ses <- options()$unmarked_shiny_session
    ses <- shiny::getDefaultReactiveDomain()
    pb <- shiny::Progress$new(ses, min=0, max=1)
    pb$set(message="Running simulations")
    if(!requireNamespace("pbapply", quietly=TRUE)){
      stop("You need to install the pbapply package", call.=FALSE)
    }
    fits <- pbapply::pblapply(1:nsim, function(i, sims, fit, bdata=NULL){
      if(!is.null(design)) fit@data <- bdata[[i]]
      if(inherits(fit, "unmarkedFitOccuMulti")){
        fit@data@ylist <- sims[[i]]
      } else{
        fit@data@y <- sims[[i]]
      }
      out <- update(fit, data=fit@data, se=TRUE)
      pb$set(value=i/nsim, message=NULL, detail=NULL)
      out
    }, sims=sims, fit=object, bdata=bdata, cl=NULL)
    pb$close()

  } else {

    fits <- lapply2(1:nsim, function(i, sims, fit, bdata=NULL){
      if(!is.null(design)) fit@data <- bdata[[i]]
      if(inherits(fit, "unmarkedFitOccuMulti")){
        fit@data@ylist <- sims[[i]]
      } else if(inherits(fit, "unmarkedFitGDR")){
        fit@data@yDistance <- sims[[i]]$yDistance
        fit@data@yRemoval <- sims[[i]]$yRemoval
      } else {
        fit@data@y <- sims[[i]]
      }
      update(fit, data=fit@data, se=TRUE)
    }, sims=sims, fit=object, bdata=bdata, cl=cl)

  }

  sum_dfs <- lapply(fits, get_summary_df)

  new("unmarkedPower", call=object@call, data=object@data, M=M,
      J=J, T=T, coefs=coefs, estimates=sum_dfs, alpha=alpha, nulls=nulls)
}

bootstrap_data <- function(data, nsims, design){
  M <- design$M
  J <- design$J
  sites <- 1:numSites(data)
  if(!is.null(J) & methods::.hasSlot(data, "numPrimary")){
    stop("Can't automatically bootstrap observations with > 1 primary period", call.=FALSE)
  }
  if(J > obsNum(data)){
    stop("Can't currently bootstrap more than the actual number of observations", call.=FALSE)
  }
  obs <- 1:obsNum(data)

  if(M > numSites(data)){
    M_samps <- lapply(1:nsims, function(i) sample(sites, M, replace=TRUE))
  } else if(M < numSites(data)){
    M_samps <- lapply(1:nsims, function(i) sample(sites, M, replace=FALSE))
  } else {
    M_samps <- replicate(nsims, sites, simplify=FALSE)
  }

  if(J > obsNum(data)){
    J_samps <- lapply(1:nsims, function(i) sample(obs, J, replace=TRUE))
  } else if(J < obsNum(data)){
    J_samps <- lapply(1:nsims, function(i) sample(obs, J, replace=FALSE))
  } else {
    J_samps <- replicate(nsims, obs, simplify=FALSE)
  }

  lapply(1:nsims, function(i) data[M_samps[[i]], J_samps[[i]]])
}

check_coefs <- function(coefs, fit, template=FALSE){
  required_subs <- names(fit@estimates@estimates)
  required_coefs <- lapply(fit@estimates@estimates, function(x) names(x@estimates))
  required_lens <- lapply(required_coefs, length)

  formulas <- sapply(names(fit), function(x) get_formula(fit, x))

  # If there are random effects, adjust the expected coefficient names
  # to remove the b vector and add the grouping covariate name
  rand <- lapply(formulas, lme4::findbars)
  if(!all(sapply(rand, is.null))){
    stopifnot(all(required_subs %in% names(formulas)))
    rvar <- lapply(rand, function(x) unlist(lapply(x, all.vars)))
    if(!all(sapply(rvar, length)<2)){
      stop("Only 1 random effect per parameter is supported", call.=FALSE)
    }
    for (i in required_subs){
      if(!is.null(rand[[i]][[1]])){
        signame <- rvar[[i]]
        old_coefs <- required_coefs[[i]]
        new_coefs <- old_coefs[!grepl("b_", old_coefs, fixed=TRUE)]
        new_coefs <- c(new_coefs, signame)
        required_coefs[[i]] <- new_coefs
      }
    }
  }

  dummy_coefs <- lapply(required_coefs, function(x){
                    out <- rep(0, length(x))
                    x <- gsub("(Intercept)", "intercept", x, fixed=TRUE)
                    names(out) <- x
                    out
                  })

  if(template) return(dummy_coefs)

  if(is.null(coefs)){
    cat("coefs argument should be a named list of named vectors, with the following structure
        (replacing 0s with your desired coefficient values):\n\n")
    print(dummy_coefs)
    stop("Supply coefs argument as specified above", call.=FALSE)
  }

  for (i in 1:length(required_subs)){
    if(!required_subs[i] %in% names(coefs)){
      stop(paste0("Missing required list element '",required_subs[i], "' in coefs list"), call.=FALSE)
    }

    sub_coefs <- coefs[[required_subs[i]]]

    if(is.null(sub_coefs)){
      stop(paste("Required coefficients for the", required_subs[i], "submodel are:",
                  paste(required_coefs[[i]],collapse=", ")))
    }

    is_named <- !is.null(names(sub_coefs)) & !any(names(sub_coefs)=="")

    if(!is_named){
      warning(paste("At least one coefficient in vector for submodel",required_subs[i],
                    "is unnamed; assuming the following order:\n",
                    paste(required_coefs[[i]], collapse=", ")))
      if(length(sub_coefs) != required_lens[i]){
        stop(paste0("Entry '",required_subs[[i]], "' in coefs list must be length ",
        required_lens[[i]]), call.=FALSE)
      }
    } else {
      rsi <- required_subs[i]
      change_int <- names(coefs[[rsi]])%in%c("intercept","Intercept")
      names(coefs[[rsi]])[change_int] <- "(Intercept)"
      change_int <- names(coefs[[rsi]])%in%c("sigmaintercept","sigmaIntercept")
      names(coefs[[rsi]])[change_int] <- "sigma(Intercept)"
      change_int <- names(coefs[[rsi]])%in%c("shapeintercept","shapeIntercept")
      names(coefs[[rsi]])[change_int] <- "shape(Intercept)"
      change_int <- names(coefs[[rsi]])%in%c("rateintercept","rateIntercept")
      names(coefs[[rsi]])[change_int] <- "rate(Intercept)"
      change_int <- grepl(" intercept", names(coefs[[rsi]]))
      names(coefs[[rsi]])[change_int] <- gsub(" intercept", " (Intercept)",
                                              names(coefs[[rsi]])[change_int])
      change_int <- grepl(" Intercept", names(coefs[[rsi]]))
      names(coefs[[rsi]])[change_int] <- gsub(" Intercept", " (Intercept)",
                                              names(coefs[[rsi]])[change_int])
      sub_coefs <- coefs[[rsi]]

      not_inc <- !required_coefs[[i]] %in% names(sub_coefs)
      extra <- !names(sub_coefs) %in% required_coefs[[i]]

      if(any(not_inc)){
        stop(paste("The following required coefficients in the", required_subs[i], "submodel were not found:",
                   paste(required_coefs[[i]][not_inc], collapse=", ")))
      }
      if(any(extra)){
        warning(paste("Ignoring extra coefficients in the", required_subs[i], "submodel:",
                      paste(names(sub_coefs)[extra], collapse=", ")))
      }
      coefs[[rsi]] <- coefs[[rsi]][required_coefs[[i]]]
    }
  }
  coefs[required_subs]
}

wald <- function(est, se, null_hyp=NULL){
  if(is.null(null_hyp) || is.na(null_hyp)) null_hyp <- 0
  Z <- (est-null_hyp)/se
  2*pnorm(abs(Z), lower.tail = FALSE)
}

diff_dir <- function(est, hyp, null_hyp=NULL){
  if(is.null(null_hyp) || is.na(null_hyp)) null_hyp <- 0
  dif <- est - null_hyp
  dif_hyp <- hyp - null_hyp
  dif * dif_hyp > 0
}

setMethod("summary", "unmarkedPower", function(object, ...){
  sum_dfs <- object@estimates
  npar <- nrow(sum_dfs[[1]])

  nulls <- object@nulls
  nulls <- lapply(nulls, function(x){
    nm <- names(x)
    nm[nm %in% c("Intercept","intercept")] <- "(Intercept)"
    names(x) <- nm
    x
  })

  coefs_no_rand <- unlist(object@coefs)[!grepl("b_", names(unlist(object@coefs)))]

  pow <- sapply(1:npar, function(ind){
    submod <- sum_dfs[[1]]$submodel[ind]
    param <- sum_dfs[[1]]$param[ind]
    ni <- nulls[[submod]][param]

    pcrit <- sapply(sum_dfs, function(x) wald(x$Estimate[ind], x$SE[ind], ni)) < object@alpha
    direct <- sapply(sum_dfs, function(x) diff_dir(x$Estimate[ind], coefs_no_rand[ind], ni))
    mean(pcrit & direct, na.rm=T)
  })

  all_nulls <- sapply(1:npar, function(ind){
    submod <- sum_dfs[[1]]$submodel[ind]
    param <- sum_dfs[[1]]$param[ind]
    ni <- nulls[[submod]][param]
    if(is.null(ni) || is.na(ni)) ni <- 0
    ni
  })

  effect_no_random <- unlist(object@coefs)[!grepl("b_",names(unlist(object@coefs)))]

  out <- cbind(sum_dfs[[1]][,1:2], effect=effect_no_random, null=all_nulls,  power=pow)
  rownames(out) <- NULL
  names(out) <- c("Submodel", "Parameter", "Effect", "Null", "Power")
  out
})

setMethod("show", "unmarkedPower", function(object){
  cat("\nModel:\n")
  print(object@call)
  cat("\n")

  cat("Power Statistics:\n")
  sumtab <- summary(object)
  sumtab$Power <- round(sumtab$Power, 3)
  print(sumtab, row.names=FALSE)
})

replace_estimates <- function(object, new_ests){
  for (i in 1:length(new_ests)){
    est <- object@estimates@estimates[[names(new_ests)[i]]]@estimates
    stopifnot(length(est) == length(new_ests[[i]]))
    object@estimates@estimates[[names(new_ests)[i]]]@estimates <- new_ests[[i]]
  }
  object
}

get_summary_df <- function(fit){
  n_est <- length(fit@estimates@estimates)
  #est_names <- unname(sapply(fit@estimates@estimates, function(x) x@name))
  est_names <- names(fit@estimates@estimates)
  all_est <- lapply(1:n_est, function(i){
    utils::capture.output(out <- summary(fit@estimates@estimates[[i]]))
    out <- cbind(submodel=est_names[i], param=rownames(out), out)
    rownames(out) <- NULL
    out
  })
  do.call(rbind, all_est)
}

setClass("unmarkedPowerList", representation(powerAnalyses="list"))

setGeneric("unmarkedPowerList", function(object, ...){
             standardGeneric("unmarkedPowerList")})

setMethod("unmarkedPowerList", "list", function(object, ...){
  new("unmarkedPowerList", powerAnalyses=object)
})

setMethod("unmarkedPowerList", "unmarkedFit",
  function(object, coefs, design, alpha=0.05, nulls=list(),
           nsim=100, parallel=FALSE, ...){

  ndesigns <- nrow(design)
  out <- lapply(1:ndesigns, function(i){
    cat(paste0("M = ",design$M[i],", J = ",design$J[i],"\n"))
    powerAnalysis(object, coefs, as.list(design[i,]), alpha=alpha, nsim=nsim,
                  nulls=nulls, parallel=FALSE)
  })
  unmarkedPowerList(out)
})

setMethod("summary", "unmarkedPowerList", function(object, ...){
  out <- lapply(object@powerAnalyses, function(x){
    stats <- summary(x)
    cbind(M=x@M, T=x@T, J=x@J, stats)
  })
  out <- do.call(rbind, out)
  out$M <- factor(out$M)
  out$T <- factor(out$T)
  out$J <- factor(out$J)
  out
})

setMethod("show", "unmarkedPowerList", function(object){
  print(summary(object))
})

setMethod("plot", "unmarkedPowerList", function(x, power=NULL, param=NULL, ...){
  dat <- summary(x)
  if(is.null(param)) param <- dat$Parameter[1]
  dat <- dat[dat$Parameter==param,,drop=FALSE]
  ylim <- range(dat$Power, na.rm=T)
  if(!is.null(power)) ylim[2] <- max(power, ylim[2])
  xlim <- range(as.numeric(as.character(dat$M)), na.rm=T)
  cols <- palette.colors(length(levels(dat$J)), palette="Dark 2")
  old_par <- graphics::par()[c("mfrow","mar")]
  nT <- length(levels(dat$T))
  mar <- old_par$mar
  if(nT == 1) mar <- c(5.1, 4.1, 2.1, 2.1)
  graphics::par(mfrow=c(length(levels(dat$T)),1), mar=mar)
  for (i in levels(dat$T)){
    plot_title <- ""
    if(nT > 1) plot_title <- paste0("T = ", i)
    tsub <- dat[dat$T==i,,drop=FALSE]
    Jlev <- levels(tsub$J)
    jsub <- tsub[tsub$J==Jlev[1],,drop=FALSE]
    plot(as.numeric(as.character(jsub$M)), jsub$Power, type="o",
        col=cols[1], ylim=ylim, xlim=xlim, xlab="Sites",
        ylab="Power", pch=19, main=plot_title)
    if(!is.null(power)) abline(h=power, lty=2)
    for (j in 2:length(Jlev)){
      jsub <- tsub[tsub$J==Jlev[j],,drop=FALSE]
      graphics::lines(as.numeric(as.character(jsub$M)), jsub$Power, type="o",
            col=cols[j], pch=19)
    }
    graphics::legend('bottomright', lwd=1, pch=19, col=cols, legend=Jlev, title="Observations")
  }
  graphics::par(mfrow=old_par)
})

setMethod("update", "unmarkedPower", function(object, ...){
  args <- list(...)
  if(!is.null(args$alpha)) object@alpha <- args$alpha
  if(!is.null(args$coefs)){
    if(!is.list(args$coefs) || all(names(args$coefs) == names(object@coefs))){
      stop("coefs list structure is incorrect", call.=FALSE)
      object@coefs <- args$coefs
    }
  }
  if(!is.null(args$nulls)) object@nulls <- args$nulls
  object
})

shinyPower <- function(object, ...){

  if(!inherits(object, "unmarkedFit")){
    stop("Requires unmarkedFit object", call.=FALSE)
  }
  if(!requireNamespace("shiny")){
    stop("Install the shiny library to use this function", call.=FALSE)
  }
  if(!requireNamespace("pbapply")){
    stop("Install the pbapply library to use this function", call.=FALSE)
  }
  options(unmarked_shiny=TRUE)
  on.exit(options(unmarked_shiny=FALSE))
  .shiny_env$.SHINY_MODEL <- object

  shiny::runApp(system.file("shinyPower", package="unmarked"))

}
