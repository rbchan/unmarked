
nmixTTD <- function(stateformula=~1, detformula=~1, data, K=100,
                       mixture=c("P","NB","ZIP"), ttdDist=c("exp", "weibull"),
                       starts, method = "BFGS",
                       se = TRUE, engine = c("C","R"), threads = 1, ...) {

  #Check arguments-------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccuTTD")){
    stop("Data is not an unmarkedFrameOccuTTD object.", call.=FALSE)
  }
  if(data@numPrimary > 1){
    stop("Multi-season data not supported.", call.=FALSE)
  }
  if(any(grepl("offset", as.character(stateformula)))){
    stop("Offsets not currently supported", call.=FALSE)
  }

  engine <- match.arg(engine)
  mixture <- match.arg(mixture)
  ttdDist <- match.arg(ttdDist)

  formula <- list(stateformula, ~1, ~1, detformula)
  check_no_support(formula)
  formula <- as.formula(paste(unlist(formula),collapse=" "))

  #Process input data----------------------------------------------------------
  designMats <- getDesign(data, formula)
  #V = detection; W = abundance
  V <- designMats$V; W <- designMats$W
  y <- designMats$y
  removed <- designMats$removed.sites

  N <- nrow(y)
  R <- ncol(y)
  T <- 1 #data@numPrimary
  J <- R / T

  #Reformat data for likelihood
  yvec <- as.numeric(t(y))
  naflag <- as.numeric(is.na(yvec))
  surveyLength <- data@surveyLength
  if(length(removed>0)) surveyLength <- surveyLength[-removed,]
  ymax <- as.numeric(t(surveyLength))
  delta <- as.numeric(yvec<ymax)

  #Organize parameters---------------------------------------------------------
  detParms <- colnames(V); nDP <- ncol(V)
  abunParms <- colnames(W); nAP <- ncol(W)

  pinds <- matrix(NA, nrow=4, ncol=2)
  pinds[1,] <- c(1, nAP)
  pinds[2,] <- c((nAP+1),(nAP+nDP))
  pinds[3,] <- nAP+nDP+1
  pinds[4,] <- nAP+nDP+(mixture%in%c("NB","ZIP"))+1

  parms <- c(abunParms, detParms)
  if(mixture == "NB") parms <- c(parms, "alpha")
  if(mixture == "ZIP") parms <- c(parms, "psi")
  if(ttdDist == "weibull") parms <- c(parms, "k")
  nP <- length(parms)

  #Likelihood functions--------------------------------------------------------

  nll_R <- function(params){

    #Get abundance and detection parameters
    lamN <- exp(W %*% params[pinds[1,]])
    lamP <- exp(V %*% params[pinds[2,]])

    if(mixture == "P"){
      pK <- sapply(0:K, function(k) dpois(k, lamN))
    } else if(mixture == "NB") {
      alpha <- exp(params[pinds[3,1]])
      pK <- sapply(0:K, function(k) dnbinom(k, mu=lamN, size = alpha))
    } else if(mixture == "ZIP"){
      alpha <- plogis(params[pinds[3,1]])
      pK <- sapply(0:K, function(k) dzip(rep(k, length(lamN)), 
                                         lambda=lamN, psi = alpha))
    }

    #Simplified version of Garrard et al. 2013 eqn 5
    #Extended to Weibull
    if(ttdDist=='weibull'){
      shape <- exp(params[pinds[4,1]])

      e_lamt <- sapply(0:K, function(k){
        lam <- k*lamP
        ( shape*lam*(lam*yvec)^(shape-1) )^delta * exp(-1*(lam*yvec)^shape)
      })

    } else {
      #Exponential
      e_lamt <- sapply(0:K, function(k) (lamP*k)^delta * exp(-lamP*k*yvec))
    }

    get_Py <- function(e_lamt, delta){
      sum_delt <- as.numeric(sum(delta, na.rm=T)>0)

      out <- rep(NA, length=K+1)
    }

    #Begin likelihood calculation
    lik <- rep(NA,N)
    ystart <- 1

    for (n in 1:N){

      yend <- ystart+J-1
      pT <- rep(NA,length=K+1)
      pT[1] <- 1 - max(delta[ystart:yend], na.rm=T)
      for (k in 1:K){
        elamt_sub <- e_lamt[ystart:yend, k+1]
        pT[k+1] <- prod(elamt_sub[!is.na(elamt_sub)])
      }
      ystart <- ystart + J

      lik[n] <- pK[n,] %*% pT
    }
    -sum(log(lik))
  }

  nll_C <- function(params){
    nll_nmixTTD(params, yvec, delta, W, V, pinds - 1, mixture, ttdDist,
                N, J, K, naflag, threads)
  }

  nll <- nll_C
  if(engine == "R") nll <- nll_R

  #Run optim()-----------------------------------------------------------------
  if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))
  if(missing(starts)) starts <- rep(0, nP)

  fm <- optim(starts, nll, method = method, hessian = se, ...)
  covMat <- invertHessian(fm, nP, se)

  #Build output object---------------------------------------------------------
  ests <- fm$par
  fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)
  names(ests) <- parms

  inds <- pinds[1,1]:pinds[1,2]
  state <- unmarkedEstimate(name = "Abundance", short.name = "lamN",
                          estimates = ests[inds],
                          covMat = as.matrix(covMat[inds,inds]),
                          invlink = "exp",
                          invlinkGrad = "exp")

  inds <- pinds[2,1]:pinds[2,2]
  det <- unmarkedEstimate(name = "Detection", short.name = "lamP",
                          estimates = ests[inds],
                          covMat = as.matrix(covMat[inds,inds]),
                          invlink = "exp",
                          invlinkGrad = "exp")


  estimateList <- unmarkedEstimateList(list(state = state, det=det))


  #Add negative binomial dispersion parameter if necessary
  if(mixture=="NB"){
    estimateList@estimates$alpha <-
      unmarkedEstimate(name = "Dispersion",
                       short.name = "alpha", estimates = ests[pinds[3,1]],
                       covMat = as.matrix(covMat[pinds[3,1], pinds[3,1]]),
                       invlink = "exp", invlinkGrad = "exp")
  }

  if(mixture=="ZIP"){
    estimateList@estimates$psi <-
      unmarkedEstimate(name = "Zero-inflation",
                       short.name = "psi", estimates = ests[pinds[3,1]],
                       covMat = as.matrix(covMat[pinds[3,1], pinds[3,1]]),
                       invlink = "logistic", invlinkGrad = "logistic.grad")
  }

  #Add Weibull shape parameter if necessary
  if(ttdDist=="weibull"){
    estimateList@estimates$shape <-
      unmarkedEstimate(name = "Weibull shape",
                       short.name = "k", estimates = ests[pinds[4,1]],
                       covMat = as.matrix(covMat[pinds[4,1], pinds[4,1]]),
                       invlink = "exp", invlinkGrad = "exp")
  }

  umfit <- new("unmarkedFitNmixTTD", fitType = "nmixTTD",
               call = match.call(),
               formula = formula,
               stateformula = stateformula,
               detformula = detformula,
               K = K,
               data = data, sitesRemoved = removed,
               estimates = estimateList,
               AIC = fmAIC, opt = fm, negLogLike = fm$value,
               nllFun = nll)

  return(umfit)
}
