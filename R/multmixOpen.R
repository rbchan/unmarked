multmixOpen <- function(lambdaformula, gammaformula, omegaformula, pformula,
    data, mixture=c("P", "NB", "ZIP"), K,
    dynamics=c("constant", "autoreg", "notrend", "trend", "ricker", "gompertz"),
    fix=c("none", "gamma", "omega"), immigration=FALSE, iotaformula = ~1,
    starts, method="BFGS", se=TRUE, ...)
{

  #Check data source
  if(!is(data, "unmarkedFrameMMO"))
    stop("Data is not of class unmarkedFrameMMO.")

  piFun <- data@piFun

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

  M <- nrow(y)
  T <- data@numPrimary
  J <- ncol(getY(data)) / T
  R <- obsNum(data) / T

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

  Xlam.offset <- D$Xlam.offset
  Xgam.offset <- D$Xgam.offset
  Xom.offset <- D$Xom.offset
  Xp.offset <- D$Xp.offset
  Xiota.offset <- D$Xiota.offset
  if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
  if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
  if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
  if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T*R)
  if(is.null(Xiota.offset)) Xiota.offset <- rep(0, M*(T-1))

  #K stuff
  if(missing(K)) {
    K <- max(y, na.rm=T) + 20
    warning("K was not specified and was set to ", K, ".")
  }
  if(K <= max(y, na.rm = TRUE))
    stop("specified K is too small. Try a value larger than any observation")
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

  lamParms <- colnames(D$Xlam)
  gamParms <- colnames(D$Xgam)
  omParms <- colnames(D$Xom)
  detParms <- colnames(D$Xp)
  nAP <- ncol(D$Xlam)
  nGP <- ncol(D$Xgam)
  nOP <- ncol(D$Xom)
  nDP <-  ncol(D$Xp)

  nIP <- ifelse(immigration, ncol(D$Xiota), 0)
  iotaParms <- character(0)
  if(immigration) iotaParms <- colnames(D$Xiota)

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

  nP <- nAP + nGP + nOP + nDP + nIP + (mixture!="P")
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))

  nbParm <- character(0)
  if(identical(mixture,"NB"))
    nbParm <- "alpha"
  else if(identical(mixture, "ZIP"))
    nbParm <- "psi"

  paramNames <- c(lamParms, gamParms, omParms, detParms,
                 iotaParms, nbParm)

  #Create indices, all possible combinations of survivors and recruits,
  #finding all unique likelihood transitions
  I <- cbind(rep(k, times=lk), rep(k, each=lk))
  I1 <- I[I[,1] <= I[,2],]
  lik_trans <- .Call("get_lik_trans", I, I1, PACKAGE="unmarked")

  beta_ind <- matrix(NA, 6, 2)
  beta_ind[1,] <- c(1, nAP) #Abundance
  beta_ind[2,] <- c(1, nGP) + nAP #Gamma
  beta_ind[3,] <- c(1, nOP) + nAP + nGP #Omega
  beta_ind[4,] <- c(1, nDP) + nAP + nGP + nOP #Detection
  beta_ind[5,] <- c(1, nIP) + nAP + nGP + nOP + nDP #Iota
  beta_ind[6,] <- c(1, 1) + nAP + nGP + nOP + nDP + nIP #2nd abun param

  #Adjustments to objects to facilitate use in c++
  fin <- fin*1 #convert to numeric
  yperm <- aperm(y, c(1,3,2))
  yna <- is.na(yperm)*1
  yna <- aperm(yna, c(3,2,1)) # fix asan problem
  yperm <- aperm(yperm, c(3,2,1)) # fix asan problem

  nll <- function(parms) {
    .Call("nll_multmixOpen",
          yperm, yt,
          D$Xlam, D$Xgam, D$Xom, D$Xp, D$Xiota,
          parms, beta_ind - 1,
          Xlam.offset, Xgam.offset, Xom.offset, Xp.offset, Xiota.offset,
          ytna, yna,
          lk, mixture, first - 1, last - 1, first1 - 1, M, T, J, R,
          D$delta, dynamics, fix, D$go.dims, immigration,
          I, I1, lik_trans$Ib, lik_trans$Ip,
          piFun, lfac.k, kmyt, lfac.kmyt, fin,
          PACKAGE = "unmarked")
  }

  if(missing(starts)){
    starts <- rep(0, nP)
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

  estimateList@estimates$det <- unmarkedEstimate(
      name = "Detection", short.name = "p",
      estimates = ests[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)],
      covMat = as.matrix(covMat[(nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP),
                        (nAP+nGP+nOP+1) : (nAP+nGP+nOP+nDP)]),
      invlink = "logistic", invlinkGrad = "logistic.grad")

  if(immigration) {
    estimateList@estimates$iota <- unmarkedEstimate(
      name="Immigration", short.name = "iota",
      estimates = ests[(nAP+nGP+nOP+nDP+1) :(nAP+nGP+nOP+nDP+nIP)],
      covMat = as.matrix(covMat[(nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP),
                                (nAP+nGP+nOP+nDP+1) : (nAP+nGP+nOP+nDP+nIP)]),
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

  umfit <- new("unmarkedFitMMO", fitType = "multmixOpen",
      call = match.call(), formula = formula, formlist = formlist, data = data,
      sitesRemoved=D$removed.sites, estimates = estimateList, AIC = fmAIC,
      opt = fm, negLogLike = fm$value, nllFun = nll, K = K, mixture = mixture,
      dynamics = dynamics, fix = fix, immigration=immigration)

  return(umfit)
}
