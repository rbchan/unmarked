setClass("unmarkedFrameGDR",
  representation(
    yDistance = "matrix",
    yRemoval = "matrix",
    survey = "character",
    dist.breaks = "numeric",
    unitsIn = "character",
    period.lengths = "numeric"
  ),
  contains="unmarkedMultFrame"
)

unmarkedFrameGDR <- function(yDistance, yRemoval, numPrimary=1,
                                     siteCovs=NULL, obsCovs=NULL,
                                     yearlySiteCovs=NULL, dist.breaks,
                                     unitsIn, period.lengths=NULL){

  if(is.null(period.lengths)){
    period.lengths <- rep(1, ncol(yRemoval)/numPrimary)
  }

  M <- nrow(yDistance)
  Jdist <- ncol(yDistance) / numPrimary
  Jrem <- ncol(yRemoval) / numPrimary

  if(length(dist.breaks) != Jdist+1){
    stop(paste("dist.breaks must have length",Jdist+1), call.=FALSE)
  }
  if(length(period.lengths) != Jrem){
    stop(paste("period.lengths must have length",Jrem), call.=FALSE)
  }

  dist_array <- array(as.vector(t(yDistance)), c(Jdist, numPrimary, M))
  dist_sums <- apply(dist_array, c(2,3), sum, na.rm=T)

  rem_array <- array(as.vector(t(yRemoval)), c(Jrem, numPrimary, M))
  rem_sums <- apply(rem_array, c(2,3), sum, na.rm=T)

  if(!all(dist_sums == rem_sums)){
    stop("Some sites/primary periods do not have the same number of distance and removal observations", call.=FALSE)
  }

  umf <- new("unmarkedFrameGDR", y=yRemoval, yDistance=yDistance,
             yRemoval=yRemoval, numPrimary=numPrimary, siteCovs=siteCovs,
             obsCovs=obsCovs, yearlySiteCovs=yearlySiteCovs, survey="point",
             dist.breaks=dist.breaks, unitsIn=unitsIn, period.lengths=period.lengths,
             obsToY=diag(ncol(yRemoval)))
  umf <- umf_to_factor(umf)
  umf
}

setAs("unmarkedFrameGDR", "data.frame", function(from){

  out <- callNextMethod(from, "data.frame")
  J <- obsNum(from)
  out <- out[,(J+1):ncol(out), drop=FALSE]

  yDistance <- from@yDistance
  colnames(yDistance) <- paste0("yDist.",1:ncol(yDistance))

  yRemoval <- from@yRemoval
  colnames(yRemoval) <- paste0("yRem.",1:ncol(yRemoval))

  data.frame(yDistance, yRemoval, out)
})

setMethod("[", c("unmarkedFrameGDR", "numeric", "missing", "missing"),
  function(x, i) {
  M <- numSites(x)
  T <- x@numPrimary

  if(length(i) == 0) return(x)
  if(any(i < 0) && any(i > 0))
    stop("i must be all positive or all negative indices.")
  if(all(i < 0)) { # if i is negative, then convert to positive
    i <- (1:M)[i]
  }

  yDist <- x@yDistance
  Rdist <- ncol(yDist)
  Jdist <- Rdist / T
  yRem <- x@yRemoval
  Rrem <- ncol(yRem)
  Jrem <- Rrem / T
  sc <- siteCovs(x)
  oc <- obsCovs(x)
  ysc <- NULL
  if(T > 1){
    ysc <- yearlySiteCovs(x)
  }

  yDist <- yDist[i,,drop=FALSE]
  yRem <- yRem[i,,drop=FALSE]

  if(!is.null(sc)){
    sc <- sc[i,,drop=FALSE]
  }

  if(!is.null(oc)){
    site_idx <- rep(1:M, each=Rrem)
    keep <- site_idx %in% i
    oc <- oc[keep,,drop=FALSE]
  }

  if(!is.null(ysc)){
    site_idx <- rep(1:M, each=T)
    keep <- site_idx %in% i
    ysc <- ysc[keep,,drop=FALSE]
  }

  umf <- x
  umf@y <- yRem
  umf@yDistance <- yDist
  umf@yRemoval <- yRem
  umf@siteCovs <- sc
  umf@obsCovs <- oc
  umf@yearlySiteCovs <- ysc

  umf
})

setMethod("[", c("unmarkedFrameGDR", "logical", "missing", "missing"),
  function(x, i) {
  i <- which(i)
  x[i, ]
})

setMethod("[", c("unmarkedFrameGDR", "missing", "numeric", "missing"),
  function(x, i, j){

  M <- numSites(x)
  T <- x@numPrimary
  if(T == 1){
    stop("Only possible to subset by primary period", call.=FALSE)
  }
  yDist <- x@yDistance
  Rdist <- ncol(yDist)
  Jdist <- Rdist / T
  yRem <- x@yRemoval
  Rrem <- ncol(yRem)
  Jrem <- Rrem / T
  oc <- obsCovs(x)
  ysc <- yearlySiteCovs(x)

  rem_idx <- rep(1:T, each=Jrem) %in% j
  yRem <- yRem[,rem_idx,drop=FALSE]
  obsToY <- x@obsToY[rem_idx, rem_idx]

  dist_idx <- rep(1:T, each=Jdist) %in% j
  yDist <- yDist[,dist_idx,drop=FALSE]

  if(!is.null(oc)){
    T_idx <- rep(rep(1:T, each=Jrem),M)
    keep <- T_idx %in% j
    oc <- oc[keep,,drop=FALSE]
  }

  if(!is.null(ysc)){
    site_idx <- rep(1:T, M)
    keep <- site_idx %in% j
    ysc <- ysc[keep,,drop=FALSE]
  }

  umf <- x
  umf@y <- yRem
  umf@yDistance <- yDist
  umf@yRemoval <- yRem
  umf@obsCovs <- oc
  umf@yearlySiteCovs <- ysc
  umf@obsToY <- obsToY
  umf@numPrimary <- length(j)

  umf
})


setMethod("getDesign", "unmarkedFrameGDR",
  function(umf, formula, na.rm=TRUE, return.frames=FALSE){

  M <- numSites(umf)
  T <- umf@numPrimary
  Rdist <- ncol(umf@yDistance)
  Jdist <- Rdist/T
  Rrem <- ncol(umf@yRemoval)
  Jrem <- Rrem/T
  yRem <- as.vector(t(umf@yRemoval))
  yDist <- as.vector(t(umf@yDistance))

  sc <- siteCovs(umf)
  oc <- obsCovs(umf)
  ysc <- yearlySiteCovs(umf)

  if(is.null(sc)) sc <- data.frame(.dummy=rep(0, M))
  if(is.null(ysc)) ysc <- data.frame(.dummy=rep(0, M*T))
  if(is.null(oc)) oc <- data.frame(.dummy=rep(0, M*Rrem))

  ysc <- cbind(ysc, sc[rep(1:M, each=T),,drop=FALSE])
  oc <- cbind(oc, ysc[rep(1:nrow(ysc), each=Jrem),,drop=FALSE])

  if(return.frames) return(list(sc=sc, ysc=ysc, oc=oc))

  lam_fixed <- remove_bars(formula$lambdaformula)
  Xlam <- model.matrix(lam_fixed,
            model.frame(lam_fixed, sc, na.action=NULL))

  phi_fixed <- remove_bars(formula$phiformula)
  Xphi <- model.matrix(phi_fixed,
            model.frame(phi_fixed, ysc, na.action=NULL))

  dist_fixed <- remove_bars(formula$distanceformula)
  Xdist <- model.matrix(dist_fixed,
            model.frame(dist_fixed, ysc, na.action=NULL))

  rem_fixed <- remove_bars(formula$removalformula)
  Xrem <- model.matrix(rem_fixed,
            model.frame(rem_fixed, oc, na.action=NULL))

  Zlam <- get_Z(formula$lambdaformula, sc)
  Zphi <- get_Z(formula$phiformula, ysc)
  Zdist <- get_Z(formula$distanceformula, ysc)
  Zrem <- get_Z(formula$removalformula, oc)
 
  # Check if there are missing yearlySiteCovs
  ydist_mat <- apply(matrix(yDist, nrow=M*T, byrow=TRUE), 1, function(x) any(is.na(x)))
  yrem_mat <- apply(matrix(yRem, nrow=M*T, byrow=TRUE), 1, function(x) any(is.na(x)))
  ok_missing_phi_covs <- ydist_mat | yrem_mat
  missing_phi_covs <- apply(Xphi, 1, function(x) any(is.na(x)))  
  if(!all(which(missing_phi_covs) %in% which(ok_missing_phi_covs))){
    stop("Missing yearlySiteCovs values for some observations that are not missing", call.=FALSE)
  }

  # Check if there are missing dist covs
  missing_dist_covs <- apply(Xdist, 1, function(x) any(is.na(x)))
  ok_missing_dist_covs <- ydist_mat
  if(!all(which(missing_dist_covs) %in% which(ok_missing_dist_covs))){
    stop("Missing yearlySiteCovs values for some distance observations that are not missing", call.=FALSE)
  }

  # Check if there are missing rem covs
  missing_obs_covs <- apply(Xrem, 1, function(x) any(is.na(x)))
  missing_obs_covs <- apply(matrix(missing_obs_covs, nrow=M*T, byrow=TRUE), 1, function(x) any(x))
  ok_missing_obs_covs <- yrem_mat
  if(!all(which(missing_obs_covs) %in% which(ok_missing_obs_covs))){
    stop("Missing obsCovs values for some removal observations that are not missing", call.=FALSE)
  }
    
  if(any(is.na(Xlam))){
    stop("gdistremoval does not currently handle missing values in siteCovs", call.=FALSE)
  }

  list(yDist=yDist, yRem=yRem, Xlam=Xlam, Xphi=Xphi, Xdist=Xdist, Xrem=Xrem,
       Zlam=Zlam, Zphi=Zphi, Zdist=Zdist, Zrem=Zrem)
})

setClass("unmarkedFitGDR", contains = "unmarkedFitGDS")

gdistremoval <- function(lambdaformula=~1, phiformula=~1, removalformula=~1,
  distanceformula=~1, data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
  output=c("abund", "density"), unitsOut=c("ha", "kmsq"), mixture=c('P', 'NB', 'ZIP'),
  K, starts, method = "BFGS", se = TRUE, engine=c("C","TMB"), threads=1, ...){

  keyfun <- match.arg(keyfun)
  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)
  mixture <- match.arg(mixture)
  engine <- match.arg(engine)

  formlist <- mget(c("lambdaformula", "phiformula", "distanceformula", "removalformula"))
  if(any(sapply(formlist, has_random))) engine <- "TMB"

  M <- numSites(data)
  T <- data@numPrimary
  Rdist <- ncol(data@yDistance)
  Rrem <- ncol(data@yRemoval)
  mixture_code <- switch(mixture, P={1}, NB={2}, ZIP={3})

  gd <- getDesign(data, formlist)

  Jdist <- Rdist / T
  ysum <- array(t(gd$yDist), c(Jdist, T, M))
  ysum <- t(apply(ysum, c(2,3), sum, na.rm=T))

  Kmin = apply(ysum, 1, max, na.rm=T)

  # Parameters-----------------------------------------------------------------
  n_param <- c(ncol(gd$Xlam), ifelse(mixture=="P",0,1),
              ifelse(T>1,ncol(gd$Xphi),0),
              ifelse(keyfun=="uniform", 0, ncol(gd$Xdist)),
              ifelse(keyfun=="hazard",1,0),
              ncol(gd$Xrem))
  nP <- sum(n_param)

  pnames <- colnames(gd$Xlam)
  if(mixture!="P") pnames <- c(pnames, "alpha")
  if(data@numPrimary > 1) pnames <- c(pnames, colnames(gd$Xphi))
  if(keyfun!="uniform") pnames <- c(pnames, colnames(gd$Xdist))
  if(keyfun=="hazard") pnames <- c(pnames, "scale")
  pnames <- c(pnames, colnames(gd$Xrem))

  lam_ind <- 1:n_param[1]
  a_ind <- n_param[1]+1
  phi_ind <- (sum(n_param[1:2])+1):(sum(n_param[1:3]))
  dist_ind <- (sum(n_param[1:3])+1):(sum(n_param[1:4]))
  sc_ind <- (sum(n_param[1:4])+1)
  rem_ind <- (sum(n_param[1:5])+1):(sum(n_param[1:6]))

  # Distance info--------------------------------------------------------------
  db <- data@dist.breaks
  w <- diff(db)
  umf_new <- data
  umf_new@y <- umf_new@yDistance
  ua <- getUA(umf_new) #in utils.R
  u <- ua$u; a <- ua$a
  A <- rowSums(a)
  switch(data@unitsIn, m = A <- A / 1e6, km = A <- A)
  switch(unitsOut,ha = A <- A * 100, kmsq = A <- A)
  if(output=='abund') A <- rep(1, numSites(data))

  # Removal info---------------------------------------------------------------
  pl <- data@period.lengths

  # Get K----------------------------------------------------------------------
  if(missing(K) || is.null(K)) K <- max(Kmin, na.rm=TRUE) + 40

  # Using C++ engine-----------------------------------------------------------
  if(engine == "C"){

    if(missing(starts)){
      starts <- rep(0, nP)
      starts[sum(n_param[1:3])+1] <- log(median(db))
    } else if(length(starts)!=nP){
      stop(paste0("starts must be length ",sum(n_param)), call.=FALSE)
    }

    nll <- function(param){
      nll_gdistremoval(param, n_param, gd$yDist, gd$yRem, ysum, mixture_code, keyfun,
                      gd$Xlam, A, gd$Xphi, gd$Xrem, gd$Xdist, db, a, t(u), w, pl,
                      K, Kmin, threads=threads)
    }

    opt <- optim(starts, nll, method=method, hessian=se, ...)

    covMat <- invertHessian(opt, nP, se)
    ests <- opt$par
    names(ests) <- pnames
    fmAIC <- 2 * opt$value + 2 * nP
    tmb_mod <- NULL

    #Organize fixed-effect estimates
    lambda_coef <- list(ests=ests[lam_ind], cov=as.matrix(covMat[lam_ind,lam_ind]))
    if(mixture != "P"){
      alpha_coef <- list(ests=ests[a_ind], cov=as.matrix(covMat[a_ind,a_ind]))
    }
    if(T > 1){
      phi_coef <- list(ests=ests[phi_ind], cov=as.matrix(covMat[phi_ind,phi_ind]))
    }
    if(keyfun != "uniform"){
      dist_coef <- list(ests=ests[dist_ind], cov=as.matrix(covMat[dist_ind,dist_ind]))
    }
    if(keyfun == "hazard"){
      scale_coef <- list(ests=ests[sc_ind], cov=as.matrix(covMat[sc_ind,sc_ind]))
    }
    rem_coef <- list(ests=ests[rem_ind], cov=as.matrix(covMat[rem_ind,rem_ind]))

    # No random effects in C engine
    lambda_rand_info <- phi_rand_info <- dist_rand_info <- rem_rand_info <- list()

  }

  # Using TMB engine-----------------------------------------------------------

  if(engine == "TMB"){
    if(missing(starts)) starts <- NULL

    dlist <- list(lambda=siteCovs(data), phi=yearlySiteCovs(data),
                  dist=siteCovs(data), rem=obsCovs(data))
    inps <- get_ranef_inputs(formlist, dlist,
                             gd[c("Xlam","Xphi","Xdist","Xrem")],
                             gd[c("Zlam","Zphi","Zdist","Zrem")])

    keyfun_type <- switch(keyfun, uniform={0}, halfnorm={1}, exp={2},
                          hazard={3})
    tmb_dat <- c(list(y_dist=gd$yDist, y_rem=gd$yRem, y_sum=ysum,
                      mixture=mixture_code, K=K, Kmin=Kmin, T=T, keyfun_type=keyfun_type,
                      A=A, db=db, a=a, w=w, u=u, per_len=pl), inps$data)

    tmb_param <- c(inps$pars, list(beta_alpha=rep(0,0), beta_scale=rep(0,0)))

    if(is.null(starts)){
      if(keyfun != "uniform") tmb_param$beta_dist[1] <- log(median(db))
    }

    if(keyfun == "hazard") tmb_param$beta_scale <- rep(0,1)
    if(keyfun == "uniform") tmb_param$beta_dist <- rep(0,0)
    if(mixture != "P") tmb_param$beta_alpha <- rep(0,1)
    if(T == 1){
      tmb_param$beta_phi <- tmb_param$b_phi <- tmb_param$lsigma_phi <- rep(0,0)
    }
    # Fit model in TMB
    tmb_out <- fit_TMB("tmb_gdistremoval", tmb_dat, tmb_param, inps$rand_ef,
                       starts=starts, method)
    tmb_mod <- tmb_out$TMB
    opt <- tmb_out$opt
    fmAIC <- tmb_out$AIC
    nll <- tmb_mod$fn

    # Organize estimates from TMB output
    lambda_coef <- get_coef_info(tmb_out$sdr, "lambda", pnames[lam_ind], lam_ind)
    lambda_rand_info <- get_randvar_info(tmb_out$sdr, "lambda", lambdaformula, siteCovs(data))
    if(mixture != "P"){
      alpha_coef <- get_coef_info(tmb_out$sdr, "alpha", pnames[a_ind], a_ind)
    }
    if(T > 1){
      phi_coef <- get_coef_info(tmb_out$sdr, "phi", pnames[phi_ind], phi_ind)
      phi_rand_info <- get_randvar_info(tmb_out$sdr, "phi", phiformula, yearlySiteCovs(data))
    }
    if(keyfun != "uniform"){
      dist_coef <- get_coef_info(tmb_out$sdr, "dist", pnames[dist_ind], dist_ind)
      dist_rand_info <- get_randvar_info(tmb_out$sdr, "dist", distanceformula, siteCovs(data))
    }
    if(keyfun == "hazard"){
      scale_coef <- get_coef_info(tmb_out$sdr, "scale", pnames[sc_ind], sc_ind)
    }
    rem_coef <- get_coef_info(tmb_out$sdr, "rem", pnames[rem_ind], rem_ind)
    rem_rand_info <- get_randvar_info(tmb_out$sdr, "rem", removalformula, obsCovs(data))

  }

  lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = lambda_coef$ests, covMat = lambda_coef$cov, fixed = 1:length(lam_ind),
    invlink = "exp", invlinkGrad = "exp", randomVarInfo=lambda_rand_info)

  estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))

  if(mixture!="P"){
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = alpha_coef$ests, covMat = alpha_coef$cov,
        fixed = 1, invlink = "exp", invlinkGrad = "exp", randomVarInfo=list())
  }

  if(T>1){
    estimateList@estimates$phi <- unmarkedEstimate(name = "Availability",
        short.name = "phi", estimates = phi_coef$ests, covMat = phi_coef$cov,
        fixed = 1:length(phi_ind), invlink = "logistic", invlinkGrad = "logistic.grad",
        randomVarInfo=phi_rand_info)
  }

  if(keyfun!="uniform"){
    estimateList@estimates$dist <- unmarkedEstimate(name = "Distance",
        short.name = "dist", estimates = dist_coef$ests, covMat = dist_coef$cov,
        fixed = 1:length(dist_ind), invlink = "exp", invlinkGrad = "exp",
        randomVarInfo=dist_rand_info)
  }

  if(keyfun=="hazard"){
    estimateList@estimates$scale <- unmarkedEstimate(name = "Hazard-rate (scale)",
        short.name = "scale", estimates = scale_coef$ests, covMat = scale_coef$cov,
        fixed = 1, invlink = "exp", invlinkGrad = "exp", randomVarInfo=list())
  }

  estimateList@estimates$rem <- unmarkedEstimate(name = "Removal",
      short.name = "rem", estimates = rem_coef$ests, covMat = rem_coef$cov,
      fixed = 1:length(rem_ind), invlink = "logistic", invlinkGrad = "logistic.grad",
      randomVarInfo=rem_rand_info)

  new("unmarkedFitGDR", fitType = "gdistremoval",
    call = match.call(), formula = as.formula(paste(formlist, collapse="")),
    formlist = formlist, data = data, estimates = estimateList, sitesRemoved = numeric(0),
    AIC = fmAIC, opt = opt, negLogLike = opt$value, nllFun = nll,
    mixture=mixture, K=K, keyfun=keyfun, unitsOut=unitsOut, output=output, TMB=tmb_mod)

}

# Methods

setMethod("getP", "unmarkedFitGDR", function(object){

  M <- numSites(object@data)
  T <- object@data@numPrimary
  Jrem <- ncol(object@data@yRemoval)/T
  Jdist <- ncol(object@data@yDistance)/T

  rem <- predict(object, "rem", level=NULL)$Predicted
  rem <- array(rem, c(Jrem, T, M))
  rem <- aperm(rem, c(3,1,2))

  pif <- array(NA, dim(rem))
  int_times <- object@data@period.lengths
  removalPiFun2 <- makeRemPiFun(int_times)
  for (t in 1:T){
    pif[,,t] <- removalPiFun2(rem[,,t])
  }

  phi <- rep(1, M*T)
  if(T>1) phi <- predict(object, "phi", level=NULL)$Predicted
  phi <- matrix(phi, M, T, byrow=TRUE)

  keyfun <- object@keyfun
  sig <- predict(object, "dist", level=NULL)$Predicted
  sig <- matrix(sig, M, T, byrow=TRUE)
  if(keyfun=="hazard") scale <- exp(coef(object, type="scale"))

  db <- object@data@dist.breaks
  a <- u <- rep(NA, Jdist)
  a[1] <- pi*db[2]^2
  for (j in 2:Jdist){
    a[j] <- pi*db[j+1]^2 - sum(a[1:(j-1)])
  }
  u <- a/sum(a)

  cp <- array(NA, c(M, Jdist, T))
  kf <- switch(keyfun, halfnorm=grhn, exp=grexp, hazard=grhaz,
               uniform=NULL)

  for (m in 1:M){
    for (t in 1:T){
      if(object@keyfun == "uniform"){
        cp[m,,t] <- u
      } else {
        for (j in 1:Jdist){
          cl <- call("integrate", f=kf, lower=db[j], upper=db[j+1], sigma=sig[m])
          names(cl)[5] <- switch(keyfun, halfnorm="sigma", exp="rate",
                                 hazard="shape")
          if(keyfun=="hazard") cl$scale=scale
          cp[m,j,t] <- eval(cl)$value * 2*pi / a[j] * u[j]
        }
      }
    }
  }

  #p_rem <- apply(pif, c(1,3), sum)
  #p_dist <- apply(cp, c(1,3), sum)

  out <- list(dist=cp, rem=pif)
  if(T > 1) out$phi <- phi
  out
})

setMethod("fitted", "unmarkedFitGDR", function(object){

  T <- object@data@numPrimary

  # Adjust log lambda when there is a random intercept
  #loglam <- log(predict(object, "lambda", level=NULL)$Predicted)
  #loglam <- E_loglam(loglam, object, "lambda")
  #lam <- exp(loglam)
  lam <- predict(object, "lambda", level=NULL)$Predicted
  if(object@output == "density"){
    ua <- getUA(object@data)
    A <- rowSums(ua$a)
    switch(object@data@unitsIn, m = A <- A / 1e6, km = A <- A)
    switch(object@unitsOut,ha = A <- A * 100, kmsq = A <- A)
    lam <- lam * A
  }

  gp <- getP(object)
  rem <- gp$rem
  dist <- gp$dist
  if(T > 1) phi <- gp$phi
  p_rem <- apply(rem, c(1,3), sum)
  p_dist <- apply(dist, c(1,3), sum)

  for (t in 1:T){
    rem[,,t] <- rem[,,t] * p_dist[,rep(t, ncol(rem[,,t]))]
    dist[,,t] <- dist[,,t] * p_rem[,rep(t,ncol(dist[,,t]))]
    if(T > 1){
      rem[,,t] <- rem[,,t] * phi[,rep(t, ncol(rem[,,t]))]
      dist[,,t] <- dist[,,t] * phi[,rep(t, ncol(dist[,,t]))]
    }
  }

  if(T > 1){
    rem_final <- rem[,,1]
    dist_final <- dist[,,1]
    for (t in 2:T){
      rem_final <- cbind(rem_final, rem[,,t])
      dist_final <- cbind(dist_final, dist[,,t])
    }
  } else {
    rem_final <- drop(rem)
    dist_final <- drop(dist)
  }

  ft_rem <- lam * rem_final
  ft_dist <- lam * dist_final
  list(dist=ft_dist, rem=ft_rem)
})

setMethod("residuals", "unmarkedFitGDR", function(object){
  ft <- fitted(object)
  list(dist=object@data@yDistance - ft$dist, rem=object@data@yRemoval-ft$rem)
})

# ranef

setMethod("ranef", "unmarkedFitGDR", function(object){

  M <- numSites(object@data)
  T <- object@data@numPrimary
  K <- object@K
  mixture <- object@mixture

  Jdist <- ncol(object@data@yDistance) / T
  ysum <- array(t(object@data@yDistance), c(Jdist, T, M))
  dist_has_na <- t(apply(ysum, c(2,3), function(x) any(is.na(x))))
  ysum <- t(apply(ysum, c(2,3), sum, na.rm=T))

  Jrem <- ncol(object@data@yRemoval) / T
  ysum_rem <- array(t(object@data@yRemoval), c(Jrem, T, M))
  rem_has_na <- t(apply(ysum_rem, c(2,3), function(x) any(is.na(x))))
  has_na <- dist_has_na | rem_has_na

  Kmin = apply(ysum, 1, max, na.rm=T)

  #loglam <- log(predict(object, "lambda", level=NULL)$Predicted)
  #loglam <- E_loglam(loglam, object, "lambda")
  #lam <- exp(loglam)
  lam <- predict(object, "lambda", level=NULL)$Predicted
  if(object@output == "density"){
    ua <- getUA(object@data)
    A <- rowSums(ua$a)
    switch(object@data@unitsIn, m = A <- A / 1e6, km = A <- A)
    switch(object@unitsOut,ha = A <- A * 100, kmsq = A <- A)
    lam <- lam * A
  }

  if(object@mixture != "P"){
    alpha <- backTransform(object, "alpha")@estimate
  }

  dets <- getP(object)
  phi <- matrix(1, M, T)
  if(T > 1){
    phi <- dets$phi
  }
  cp <- dets$dist
  pif <- dets$rem

  pr <- apply(cp, c(1,3), sum)
  prRem <- apply(pif, c(1,3), sum)

  post <- array(0, c(M, K+1, 1))
  colnames(post) <- 0:K
  for (i in 1:M){
    if(mixture=="P"){
      f <- dpois(0:K, lam[i])
    } else if(mixture=="NB"){
      f <- dnbinom(0:K, mu=lam[i], size=alpha)
    } else if(mixture=="ZIP"){
      f <- dzip(0:K, lam[i], alpha)
    }

    # All sampling periods at site i have at least one missing value
    if(all(has_na[i,])){
      g <- rep(NA,K+1)
      next
    } else {
      # At least one sampling period wasn't missing
      g <- rep(1, K+1)
      for (t in 1:T){
        if(has_na[i,t]){
          next
        }
        for (k in 1:(K+1)){
          g[k] <- g[k] * dbinom(ysum[i,t], k-1, prob=pr[i,t]*prRem[i,t]*phi[i,t],
                                log=FALSE)
        }
      }
    }
    fg <- f*g
    post[i,,1] <- fg/sum(fg)
  }

  new("unmarkedRanef", post=post)
})

setMethod("simulate", "unmarkedFitGDR", function(object, nsim, seed=NULL, na.rm=FALSE){

  # Adjust log lambda when there is a random intercept
  #loglam <- log(predict(object, "lambda", level=NULL)$Predicted)
  #loglam <- E_loglam(loglam, object, "lambda")
  #lam <- exp(loglam)
  lam <- predict(object, "lambda", level=NULL)$Predicted
  if(object@output == "density"){
    ua <- getUA(object@data)
    A <- rowSums(ua$a)
    switch(object@data@unitsIn, m = A <- A / 1e6, km = A <- A)
    switch(object@unitsOut,ha = A <- A * 100, kmsq = A <- A)
    lam <- lam * A
  }
  dets <- getP(object)

  if(object@mixture != "P"){
    alpha <- backTransform(object, "alpha")@estimate
  }

  M <- length(lam)
  T <- object@data@numPrimary

  if(T > 1){
    phi <- dets$phi
  } else {
    phi <- matrix(1, M, T)
  }

  Jrem <- dim(dets$rem)[2]
  Jdist <- dim(dets$dist)[2]

  p_dist <- apply(dets$dist, c(1,3), sum)
  p_rem <- apply(dets$rem, c(1,3), sum)

  dist_scaled <- array(NA, dim(dets$dist))
  rem_scaled <- array(NA, dim(dets$rem))
  for (t in 1:T){
    dist_scaled[,,t] <- dets$dist[,,t] / p_dist[,t]
    rem_scaled[,,t] <- dets$rem[,,t] / p_rem[,t]
  }

  p_total <- p_dist * p_rem * phi
  stopifnot(dim(p_total) == c(M, T))

  out <- vector("list", nsim)

  for (i in 1:nsim){

    switch(object@mixture,
      P = N <- rpois(M, lam),
      NB = N <- rnbinom(M, size=alpha, mu=lam),
      ZIP = N <- rzip(M, lam, alpha)
    )

    ydist <- matrix(NA, M, T*Jdist)
    yrem <- matrix(NA, M, T*Jrem)

    for (m in 1:M){
      ysum <- suppressWarnings(rbinom(T, N[m], p_total[m,]))

      ydist_m <- yrem_m <- c()

      for (t in 1:T){
        if(is.na(ysum[t])){
          yrem_m <- c(yrem_m, rep(NA, Jrem))
          ydist_m <- c(ydist_m, rep(NA, Jdist))
        } else {
          rem_class <- sample(1:Jrem, ysum[t], replace=TRUE, prob=rem_scaled[m,,t])
          rem_class <- factor(rem_class, levels=1:Jrem)
          yrem_m <- c(yrem_m, as.numeric(table(rem_class)))
          dist_class <- sample(1:Jdist, ysum[t], replace=TRUE, prob=dist_scaled[m,,t])
          dist_class <- factor(dist_class, levels=1:Jdist)
          ydist_m <- c(ydist_m, as.numeric(table(dist_class)))
        }
      }
      stopifnot(length(ydist_m)==ncol(ydist))
      stopifnot(length(yrem_m)==ncol(yrem))

      ydist[m,] <- ydist_m
      yrem[m,] <- yrem_m
    }
    out[[i]] <- list(yRemoval=yrem, yDistance=ydist)
  }
  out
})


setMethod("update", "unmarkedFitGDR",
    function(object, lambdaformula, phiformula, removalformula, distanceformula,
             ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")

    if(!missing(lambdaformula)){
      call$lambdaformula <- lambdaformula
    }
    if(!missing(phiformula)){
      if(!is.null(call$phiformula)){
        call$phiformula <- phiformula
      }
    }
    if(!missing(removalformula)){
      call$removalformula <- removalformula
    }
    if(!missing(distanceformula)){
      call$distanceformula <- distanceformula
    }

    extras <- match.call(call=sys.call(-1),
                         expand.dots = FALSE)$...
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing])
            call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
            }
        }
    if (evaluate)
        eval(call, parent.frame(2))
    else call
})


setMethod("replaceY", "unmarkedFrameGDR",
          function(object, newY, replNA=TRUE, ...){

      ydist <- newY$yDistance
      stopifnot(dim(ydist)==dim(object@yDistance))
      yrem <- newY$yRemoval
      stopifnot(dim(yrem)==dim(object@yRemoval))

      if(replNA){
        ydist[is.na(object@yDistance)] <- NA
        yrem[is.na(object@yRemoval)] <- NA
      }

      object@yDistance <- ydist
      object@yRemoval <- yrem
      object
})


setMethod("SSE", "unmarkedFitGDR", function(fit, ...){
    r <- sapply(residuals(fit), function(x) sum(x^2, na.rm=T))
    return(c(SSE = sum(r)))
})

setMethod("nonparboot", "unmarkedFitGDR",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
   stop("Not currently supported for unmarkedFitGDR", call.=FALSE)
})


setMethod("plot", c(x = "unmarkedFitGDR", y = "missing"), function(x, y, ...)
{
    r <- residuals(x)
    e <- fitted(x)

    old_mfrow <- graphics::par("mfrow")
    on.exit(graphics::par(mfrow=old_mfrow))
    graphics::par(mfrow=c(2,1))

    plot(e[[1]], r[[1]], ylab="Residuals", xlab="Predicted values",
         main="Distance")
    abline(h = 0, lty = 3, col = "gray")

    plot(e[[2]], r[[2]], ylab="Residuals", xlab="Predicted values",
         main="Removal")
    abline(h = 0, lty = 3, col = "gray")
})

# Used with fitList
setMethod("fl_getY", "unmarkedFitGDR", function(fit, ...){
  getDesign(getData(fit), fit@formlist)$yDist
})
