setClass("unmarkedFitGOccu",
    representation(
        formlist = "list"),
    contains = "unmarkedFit")

setClass("unmarkedFrameGOccu", contains = "unmarkedFrameG3")

setMethod("getDesign", "unmarkedFrameGOccu",
  function(umf, formula, na.rm=TRUE){
  out <- methods::callNextMethod(umf, formula=formula, na.rm=na.rm)
  names(out)[2] <- "Xpsi"
  names(out)[5] <- "Xpsi.offset"
  out
})

unmarkedFrameGOccu <- function(y, siteCovs=NULL, obsCovs=NULL, numPrimary,
                             yearlySiteCovs=NULL) {
  y[y > 1] <- 1
  if(numPrimary < 2) stop("numPrimary < 2, use occu instead")
  umf <- unmarkedFrameGPC(y, siteCovs=siteCovs, obsCovs=obsCovs, 
                          numPrimary=numPrimary, yearlySiteCovs=yearlySiteCovs)
  class(umf) <- "unmarkedFrameGOccu"
  umf
}

goccu <- function(psiformula, phiformula, pformula, data,
                  linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                  se = TRUE, ...){

  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))
  psiLinkFunc <- ifelse(linkPsi=="cloglog", cloglog, plogis)
  psiInvLink <- ifelse(linkPsi=="cloglog", "cloglog", "logistic")
  psiLinkGrad <- ifelse(linkPsi=="cloglog", "cloglog.grad", "logistic.grad")

  # Pass phiformula as gamma/eps formula so it will be applied to
  # yearlySiteCovs in getDesign
  formlist <- list(psiformula=psiformula, phi=phiformula,
                  pformula=pformula)

  formula <- as.formula(paste(unlist(formlist), collapse=" "))

  data@y[data@y > 1] <- 1
 
  class(data) <- "unmarkedFrameGOccu"

  # handle offsets

  gd <- getDesign(data, formula = formula)

  y <- gd$y
  Xpsi <- gd$Xpsi
  Xphi <- gd$Xphi
  Xp <- gd$Xdet

  M <- nrow(y)
  T <- data@numPrimary
  J <- ncol(y) / T

  # Identify entirely missing primary periods at each site
  y_array <- array(t(y), c(J, T, M)) 
  missing_session <- t(apply(y_array, c(2,3), 
                           function(x) as.numeric(all(is.na(x)))))

  # Create possible states in each T
  alpha_potential <- as.matrix(expand.grid(rep(list(c(0, 1)), T)))
  n_possible <- nrow(alpha_potential)

  # Known present at each site
  known_present <- rep(0, M)
  # Known available at each site and T
  known_available <- matrix(0, nrow=M, ncol=T)
  
  for (i in 1:M){
    for (t in 1:T){
      for (j in 1:J){
        if(is.na(y_array[j,t,i])) next
        if(y_array[j, t, i] == 1){
          known_present[i] <- 1
          known_available[i,t] <- 1
        }
      }
    }
  }

  # Bundle data for TMB
  dataList <- list(y=y, T=T, link=ifelse(linkPsi=='cloglog', 1, 0), 
                   Xpsi=Xpsi, Xphi=Xphi, Xp=Xp,
                   n_possible=n_possible,
                   alpha_potential=alpha_potential,
                   known_present=known_present, known_available=known_available, 
                   missing_session=missing_session)

  # Provide dimensions and starting values for parameters
  # This part should change to be more like occu() if we add random effects
  psi_ind <- 1:ncol(Xpsi)
  phi_ind <- 1:ncol(Xphi) + max(psi_ind)
  p_ind <- 1:ncol(Xp) + max(phi_ind)
  nP <- max(p_ind)
  params <- list(beta_psi = rep(0, length(psi_ind)), 
                 beta_phi = rep(0, length(phi_ind)), 
                 beta_p = rep(0, length(p_ind)))

  # Create TMB object
  tmb_mod <- TMB::MakeADFun(data = c(model = "tmb_goccu", dataList),
                            parameters = params,
                            DLL = "unmarked_TMBExports", silent = TRUE)

  # Optimize TMB object, print and save results
  if(missing(starts) || is.null(starts)) starts <- tmb_mod$par
  opt <- optim(starts, fn = tmb_mod$fn, gr = tmb_mod$gr, method = method,
               hessian = se, ...)

  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  names(ests) <- c(colnames(Xpsi), colnames(Xphi), colnames(Xp))
  fmAIC <- 2 * opt$value + 2 * nP


  psi_est <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                              estimates = ests[psi_ind],
                              covMat = covMat[psi_ind, psi_ind, drop=FALSE],
                              fixed = 1:ncol(Xpsi),
                              invlink = psiInvLink,
                              invlinkGrad = psiLinkGrad,
                              randomVarInfo=list()
                            )

  phi_est <- unmarkedEstimate(name = "Availability", short.name = "phi",
                              estimates = ests[phi_ind],
                              covMat = covMat[phi_ind, phi_ind, drop=FALSE],
                              fixed = 1:ncol(Xphi),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad",
                              randomVarInfo=list()
                            )

  p_est <- unmarkedEstimate(name = "Detection", short.name = "p",
                              estimates = ests[p_ind],
                              covMat = covMat[p_ind, p_ind, drop=FALSE],
                              fixed = 1:ncol(Xp),
                              invlink = "logistic",
                              invlinkGrad = "logistic.grad",
                              randomVarInfo=list()
                            )

  estimate_list <- unmarkedEstimateList(list(psi=psi_est, phi=phi_est,
                                                        det=p_est))

  # Create unmarkedFit object--------------------------------------------------
  umfit <- new("unmarkedFitGOccu", fitType = "goccu", call = match.call(),
                 formula = formula, formlist=formlist, data = data,
                 sitesRemoved = gd$removed.sites,
                 estimates = estimate_list, AIC = fmAIC, opt = opt,
                 negLogLike = opt$value,
                 nllFun = tmb_mod$fn, TMB=tmb_mod)

  return(umfit)

}

# Methods

setMethod("predict_inputs_from_umf", "unmarkedFitGOccu",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm=na.rm)
  X_idx <- switch(type, psi="Xpsi", phi="Xphi", det="Xdet")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitGOccu", function(object, type, ...){
  fl <- object@formlist
  switch(type, psi=fl$psiformula, phi=fl$phiformula, det=fl$pformula)
})

setMethod("get_orig_data", "unmarkedFitGOccu", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=FALSE)
  datatype <- switch(type, psi='site_covs', phi='yearly_site_covs',
                     det='obs_covs')
  clean_covs[[datatype]]
})

setMethod("getP", "unmarkedFitGOccu",
  function(object, na.rm=FALSE){
  gd <- getDesign(object@data, object@formula, na.rm=na.rm)
  p <- drop(plogis(gd$Xdet %*% coef(object, "det")))
  M <- numSites(object@data)
  p <- matrix(p, nrow=M, ncol=obsNum(object@data), 
              byrow=TRUE)
  p
})

setMethod("fitted", "unmarkedFitGOccu",
  function(object, na.rm= FALSE){

  M <- numSites(object@data)
  JT <- obsNum(object@data)  
  gd <- getDesign(object@data, object@formula, na.rm=na.rm)

  psi <- drop(plogis(gd$Xpsi %*% coef(object, "psi")))
  psi <- matrix(psi, nrow=M, ncol=JT)

  phi <- drop(plogis(gd$Xphi %*% coef(object, "phi")))
  phi <- rep(phi, each = JT / object@data@numPrimary)
  phi <- matrix(phi, nrow=M, ncol=JT, byrow=TRUE)

  p <- getP(object)

  psi * phi * p
})


# based on ranef for GPC
setMethod("ranef", "unmarkedFitGOccu", function(object, ...){

  M <- numSites(object@data)
  JT <- obsNum(object@data)
  T <- object@data@numPrimary
  J <- JT / T

  gd <- getDesign(object@data, object@formula, na.rm=FALSE)
  y_array <- array(t(gd$y), c(J, T, M)) 

  psi <- drop(plogis(gd$Xpsi %*% coef(object, "psi")))
  phi <- drop(plogis(gd$Xphi %*% coef(object, "phi")))
  phi <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
  p <- getP(object)
  p_array <- array(t(p), c(J, T, M))
  
  Z <- ZZ <- 0:1
  post <- array(0, c(M, 2, 1))
  colnames(post) <- Z

  for(i in 1:M) {
    f <- dbinom(Z, 1, psi[i])
    
    ghi <- rep(0, 2)

    for(t in 1:T) {
      gh <- matrix(-Inf, 2, 2)
      for(z in Z) {
        if(z < max(y_array[,,i], na.rm=TRUE)){
          gh[,z+1] <- -Inf
          next
        }
        if(is.na(phi[i,t])) {
          g <- rep(0, 2)
          g[ZZ>z] <- -Inf
        } else{
          g <- dbinom(ZZ, z, phi[i,t], log=TRUE)
        }
        h <- rep(0, 2)
        for(j in 1:J) {
          if(is.na(y_array[j,t,i]) | is.na(p_array[j,t,i])) next
          h <- h + dbinom(y_array[j,t,i], ZZ, p_array[j,t,i], log=TRUE)
        }
        gh[,z+1] <- g + h
      }
      ghi <- ghi + log(colSums(exp(gh)))
    }
    fgh <- exp(f + ghi)
    prM <- fgh/sum(fgh)
    post[i,,1] <- prM
  }

  new("unmarkedRanef", post=post)
})


setMethod("simulate", "unmarkedFitGOccu", 
          function(object, nsim = 1, seed = NULL, na.rm = FALSE){
  
  gd <- getDesign(object@data, object@formula, na.rm=FALSE)
  M <- nrow(gd$y)
  T <- object@data@numPrimary
  JT <- ncol(gd$y)
  J <- JT / T
  y_array <- array(t(gd$y), c(J, T, M)) 

  psi <- drop(plogis(gd$Xpsi %*% coef(object, "psi")))
  phi <- drop(plogis(gd$Xphi %*% coef(object, "phi")))
  phi <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
  p <- getP(object)

  sim_list <- list()

  for (i in 1:nsim){
    z <- suppressWarnings(rbinom(M, 1, psi))
    z <- matrix(z, nrow=M, ncol=T) 
    
    zz <- suppressWarnings(rbinom(M*T, 1, phi*z))
    zz <- matrix(zz, M, T)
    
    colrep <- rep(1:T, each=J)
    zz <- zz[,colrep]

    y <- suppressWarnings(rbinom(M*T*J, 1, zz*p))
    y <- matrix(y, M, JT)
    if(na.rm) y[which(is.na(gd$y))] <- NA 
    sim_list[[i]] <- y
  }

  return(sim_list)
})


setMethod("update", "unmarkedFitGOccu",
    function(object, psiformula, phiformula, pformula, ...,
        evaluate = TRUE)
{
    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    formlist <- object@formlist
    if (!missing(psiformula))
        call$psiformula <- update.formula(formlist$psiformula, psiformula)
    if (!missing(phiformula))
        call$phiformula <- update.formula(formlist$phiformula, phiformula)
    if (!missing(pformula))
        call$pformula <- update.formula(formlist$pformula, pformula)
    extras <- match.call(call=sys.call(-1),
                         expand.dots = FALSE)$...
    if(length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
            }
        }
    if (evaluate)
        eval(call, parent.frame(2))
    else call
})
