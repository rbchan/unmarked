distsampOpen <- function(lambdaformula, gammaformula, omegaformula, pformula,
    data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("abund", "density"), unitsOut=c("ha", "kmsq"),
    mixture=c("P", "NB", "ZIP"), K,
    dynamics=c("constant", "autoreg", "notrend", "trend", "ricker", "gompertz"),
    fix=c("none", "gamma", "omega"), immigration=FALSE, iotaformula = ~1,
    starts, method="BFGS", se=TRUE, ...)
{

  #Check data source
  if(!is(data, "unmarkedFrameDSO"))
    stop("Data is not of class unmarkedFrameDSO.")

  #Check detection model arguments
  keyfun <- match.arg(keyfun)
  if(!keyfun %in% c("halfnorm", "exp", "hazard", "uniform"))
    stop("keyfun must be 'halfnorm', 'exp', 'hazard', or 'uniform'")
  if(keyfun == "uniform"){
    if(!missing(pformula)){
      warning("pformula is ignored when using a uniform key function")
    }
    pformula <- ~1
  }

  output <- match.arg(output)
  unitsOut <- match.arg(unitsOut)

  db <- data@dist.breaks
  w <- diff(db)
  tlength <- data@tlength
  survey <- data@survey
  unitsIn <- data@unitsIn

  #Check state model arguments
  mixture <- match.arg(mixture)
  dynamics <- match.arg(dynamics)

  if((identical(dynamics, "constant") || identical(dynamics, "notrend")) & immigration)
    stop("You can not include immigration in the constant or notrend models")

  if(identical(dynamics, "notrend") &
   !identical(lambdaformula, omegaformula))
    stop("lambdaformula and omegaformula must be identical for notrend model")

  fix <- match.arg(fix)

  formlist <- mget(c("lambdaformula", "gammaformula", "omegaformula",
                   "pformula", "iotaformula"))
  check_no_support(formlist)
  formula <- as.formula(paste(unlist(formlist), collapse=" "))

  D <- getDesign(data, formula)
  y <- D$y

  Xlam <- D$Xlam
  Xgam <- D$Xgam
  Xom <- D$Xom
  Xsig <- D$Xp
  Xiota<- D$Xiota

  delta <- D$delta; go.dims <- D$go.dims
  deltamax <- max(delta, na.rm=TRUE)
  M <- nrow(y)
  T <- data@numPrimary
  J <- ncol(getY(data)) / T

  Xlam.offset <- D$Xlam.offset
  Xgam.offset <- D$Xgam.offset
  Xom.offset <- D$Xom.offset
  Xsig.offset <- D$Xp.offset
  Xiota.offset<- D$Xiota.offset

  y <- array(y, c(M, J, T))
  yt <- apply(y, c(1,3), function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
  })

  ytna <- apply(is.na(y), c(1,3), all)
  ytna <- matrix(ytna, nrow=M)
  ytna[] <- as.integer(ytna)

  first <- apply(!ytna, 1, function(x) min(which(x)))
  last  <- apply(!ytna, 1, function(x) max(which(x)))
  first1 <- which(first==1)[1]

  #K stuff
  if(missing(K)) {
    K <- max(y, na.rm=T) + 20
    warning("K was not specified and was set to ", K, ".")
  }

  J <- ncol(data@y) / data@numPrimary
  inds <- split(1:ncol(data@y), ceiling(1:ncol(data@y)/J))
  Tobs <- sapply(1:length(inds), function(i){
    rowSums(data@y[,inds[[i]], drop=FALSE], na.rm=TRUE)
  })
  Kmin <- max(Tobs, na.rm=TRUE)

  if(K < Kmin){
    stop("Specified K is too small, must be larger than the max total count in a primary period",
         call.=FALSE)
  }
  k <- 0:K
  lk <- length(k)
  #Some k-related indices to avoid repeated calculations in likelihood
  lfac.k <- lgamma(k+1)
  kmyt <- array(0, c(lk, T, M))
  lfac.kmyt <- array(0, c(M, T, lk))
  fin <- array(NA, c(M, T, lk)) #Indicator if given k is possible given y
  for(i in 1:M) {
    for(t in 1:T) {
      fin[i,t,] <- k - yt[i,t] >= 0
      if(sum(ytna[i,t])==0) {
        kmyt[,t,i] <- k - yt[i,t]
        lfac.kmyt[i,t, ] <- lgamma(kmyt[,t,i] + 1)
      }
    }
  }

  #Transect areas / proportions
  ua <- getUA(data) #in utils.R
  u <- ua$u; a <- ua$a
  if(length(D$removed.sites)>0){
    u <- ua$u[-D$removed.sites,]
    a <- ua$a[-D$removed.sites,]
  }

  switch(survey,
    line = A <- rowSums(a) * 2,
    point = A <- rowSums(a))
  switch(unitsIn,
    m = A <- A / 1e6,
    km = A <- A)
  switch(unitsOut,
    ha = A <- A * 100,
    kmsq = A <- A)
  if(output=='abund'){
    A <- rep(1, M)
  }

  lamParms <- colnames(Xlam)
  gamParms <- colnames(Xgam)
  omParms <- colnames(Xom)
  nAP <- ncol(Xlam)
  nGP <- ncol(Xgam)
  nOP <- ncol(Xom)

  #No parameters if uniform key function
  nDP <- ifelse(keyfun == "uniform", 0, ncol(Xsig))
  detParms <- character(0)
  if(keyfun != "uniform") detParms <- colnames(Xsig)

  nIP <- ifelse(immigration, ncol(Xiota), 0)
  iotaParms <- character(0)
  if(immigration) iotaParms <- colnames(Xiota)

  if(identical(fix, "gamma")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nGP > 1){
        stop("gamma covariates not allowed when fix==gamma")
    }else {
        nGP <- 0
        gamParms <- character(0)
    }
  } else if(identical(dynamics, "notrend")) {
    if(nGP > 1){
        stop("gamma covariates not allowed when dyamics==notrend")
    } else {
        nGP <- 0
        gamParms <- character(0)
    }
  }

  if(identical(fix, "omega")) {
    if(!identical(dynamics, "constant"))
        stop("dynamics must be constant when fixing gamma or omega")
    if(nOP > 1)
        stop("omega covariates not allowed when fix==omega")
    else {
        nOP <- 0
        omParms <- character(0)
    }
  } else if(identical(dynamics, "trend")) {
    if(nOP > 1)
        stop("omega covariates not allowed when dynamics='trend'")
    else {
        nOP <- 0
        omParms <- character(0)
    }
  }

  nP <- nAP + nGP + nOP + nDP + nIP + (mixture!="P") + (keyfun == "hazard")
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))

  nbParm <- character(0)
  if(identical(mixture,"NB"))
    nbParm <- "alpha"
  else if(identical(mixture, "ZIP"))
    nbParm <- "psi"

  scaleParm <- character(0)
  if(identical(keyfun, "hazard")) scaleParm <- "scale"

  paramNames <- c(lamParms, gamParms, omParms, detParms,
                 iotaParms, scaleParm, nbParm)

  #Create indices, all possible combinations of survivors and recruits,
  #finding all unique likelihood transitions
  I <- cbind(rep(k, times=lk), rep(k, each=lk))
  I1 <- I[I[,1] <= I[,2],]
  lik_trans <- get_lik_trans(I, I1)

  beta_ind <- matrix(NA, 7, 2)
  beta_ind[1,] <- c(1, nAP) #Abundance
  beta_ind[2,] <- c(1, nGP) + nAP #Gamma
  beta_ind[3,] <- c(1, nOP) + nAP + nGP #Omega
  beta_ind[4,] <- c(1, nDP) + nAP + nGP + nOP #Sigma
  beta_ind[5,] <- c(1, nIP) + nAP + nGP + nOP + nDP #Iota
  beta_ind[6,] <- c(1, 1) + nAP + nGP + nOP + nDP + nIP #Hazard scale
  beta_ind[7,] <- c(1, 1) + nAP + nGP + nOP + nDP + nIP + (keyfun == "hazard")

  #Adjustments to objects to facilitate use in c++
  fin <- fin*1 #convert to numeric
  u <- t(u) #easier to access column-wise
  yperm <- aperm(y, c(1,3,2))

  nll <- function(parms) {
    nll_distsampOpen(
          yperm, yt,
          Xlam, Xgam, Xom, Xsig, Xiota,
          parms, beta_ind - 1,
          Xlam.offset, Xgam.offset, Xom.offset, Xsig.offset, Xiota.offset,
          ytna,
          lk, mixture, first - 1, last - 1, first1 - 1, M, T,
          delta, dynamics, survey, fix, go.dims, immigration,
          I, I1, lik_trans$Ib, lik_trans$Ip,
          a, u, w, db,
          keyfun, lfac.k, kmyt, lfac.kmyt, fin, A
          )
  }

  if(missing(starts)){
    starts <- rep(0, nP)
    #Need a semi-realistic value for sigma intercept
    if(keyfun != "uniform")
      starts[beta_ind[4,1]] <- log(mean(w))
  }

  fm <- optim(starts, nll, method=method, hessian=se, ...)
  ests <- fm$par
  names(ests) <- paramNames
  covMat <- invertHessian(fm, nP, se)
  fmAIC <- 2*fm$value + 2*nP

  lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lam",
                    estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP,1:nAP]),
                    invlink = "exp", invlinkGrad = "exp")
  estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))

  gamName <- switch(dynamics, constant = "gamConst", autoreg = "gamAR",
                              notrend = "", trend = "gamTrend",
                              ricker="gamRicker", gompertz = "gamGomp")
  if(!(identical(fix, "gamma") | identical(dynamics, "notrend"))){
    estimateList@estimates$gamma <- unmarkedEstimate(name =
        ifelse(identical(dynamics, "constant") | identical(dynamics, "autoreg"),
        "Recruitment", "Growth Rate"), short.name = gamName,
        estimates = ests[(nAP+1) : (nAP+nGP)], covMat = as.matrix(covMat[(nAP+1) :
                           (nAP+nGP), (nAP+1) : (nAP+nGP)]),
        invlink = "exp", invlinkGrad = "exp")
  }

  if(!(identical(fix, "omega") | identical(dynamics, "trend"))) {
    if(identical(dynamics, "constant") | identical(dynamics, "autoreg") |
       identical(dynamics, "notrend")){
        estimateList@estimates$omega <- unmarkedEstimate( name="Apparent Survival",
          short.name = "omega", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
          covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
                                    (nAP+nGP+1) : (nAP+nGP+nOP)]),
          invlink = "logistic", invlinkGrad = "logistic.grad")
    } else if(identical(dynamics, "ricker")){
        estimateList@estimates$omega <- unmarkedEstimate(name="Carrying Capacity",
          short.name = "omCarCap", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
          covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
                            (nAP+nGP+1) : (nAP+nGP+nOP)]),
          invlink = "exp", invlinkGrad = "exp")
    } else{
      estimateList@estimates$omega <- unmarkedEstimate(name="Carrying Capacity",
        short.name = "omCarCap", estimates = ests[(nAP+nGP+1) :(nAP+nGP+nOP)],
        covMat = as.matrix(covMat[(nAP+nGP+1) : (nAP+nGP+nOP),
                                  (nAP+nGP+1) : (nAP+nGP+nOP)]),
        invlink = "exp", invlinkGrad = "exp")
    }
  }

  if(keyfun != "uniform"){
    estimateList@estimates$det <- unmarkedEstimate(
      name = "Detection", short.name = "sigma",
      estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)],
      covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP),
                        (nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
      invlink = "exp", invlinkGrad = "exp")
  }

  if(immigration) {
    estimateList@estimates$iota <- unmarkedEstimate(
      name="Immigration", short.name = "iota",
      estimates = ests[(nAP+nGP+nOP+nDP+1) :(nAP+nGP+nOP+nDP+nIP)],
      covMat = as.matrix(covMat[(nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP),
                                (nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP)]),
      invlink = "exp", invlinkGrad = "exp")
  }

  if(identical(keyfun, "hazard")) {
    estimateList@estimates$scale <- unmarkedEstimate(name = "Hazard-rate(scale)",
        short.name = "scale", estimates = ests[nAP+nGP+nOP+nDP+nIP+1],
        covMat = as.matrix(covMat[nAP+nGP+nOP+nDP+nIP+1,
                                  nAP+nGP+nOP+nDP+nIP+1]),
        invlink = "exp", invlinkGrad = "exp")
  }

  if(identical(mixture, "NB")) {
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")
  }
  if(identical(mixture, "ZIP")) {
    estimateList@estimates$psi <- unmarkedEstimate(name = "Zero-inflation",
        short.name = "psi", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
  }

  umfit <- new("unmarkedFitDSO", fitType = "distsampOpen",
      call = match.call(), formula = formula, formlist = formlist, data = data,
      sitesRemoved=D$removed.sites, estimates = estimateList, AIC = fmAIC,
      opt = fm, negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture,
      dynamics = dynamics, fix = fix, immigration=immigration, keyfun=keyfun,
      unitsOut=unitsOut)

  return(umfit)
}
