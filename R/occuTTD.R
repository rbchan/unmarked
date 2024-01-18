
#  Fit the time-to-detection model of Garrard et al. 2008, 2013
#  with dynamic extension

occuTTD <- function(psiformula=~1, gammaformula=~1, epsilonformula=~1,
                    detformula=~1, data,  ttdDist=c('exp','weibull'),
                    linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                    se = TRUE, engine = c("C", "R"), ...) {

  #Check arguments-------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccuTTD")){
    stop("Data is not an unmarkedFrameOccuTTD object.")
  }

  engine <- match.arg(engine, c("C", "R"))
  ttdDist <- match.arg(ttdDist, c("exp","weibull"))
  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))

  formula <- list(psiformula, gammaformula, epsilonformula, detformula)
  check_no_support(formula)
  formula <- as.formula(paste(unlist(formula),collapse=" "))

  #Psi link function
  linkFunc <- plogis
  invlink <- "logistic"
  linkGrad <- "logistic.grad"
  if(linkPsi == "cloglog"){
    linkFunc <- cloglog
    invlink <- "cloglog"
    linkGrad <- "cloglog.grad"
  }

  #Process input data----------------------------------------------------------
  designMats <- getDesign(data, formula)
  #V = detection; W = occupancy
  V <- designMats$V; W <- designMats$W
  X.gam <- designMats$X.gam; X.eps <- designMats$X.eps
  y <- designMats$y
  removed <- designMats$removed.sites

  N <- nrow(y)
  R <- ncol(y)
  T <- data@numPrimary
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
  occParms <- colnames(W); nOP <- ncol(W)
  psi_inds <- 1:nOP

  gamParms <- NULL; nGP <- 0; col_inds <- c(0,0)
  epsParms <- NULL; nEP <- 0; ext_inds <- c(0,0)
  if(T>1){

    #Remove final year from col/ext design matrices
    X.gam <- as.matrix(X.gam[-seq(T,N*T,by=T),,drop=FALSE])
    X.eps <- as.matrix(X.eps[-seq(T,N*T,by=T),,drop=FALSE])

    gamParms <- colnames(X.gam); nGP <- ncol(X.gam)
    epsParms <- colnames(X.eps); nEP <- ncol(X.eps)
    col_inds <- (nOP+1):(nOP+nGP)
    ext_inds <- (nOP+nGP+1):(nOP+nGP+nEP)
  }
  det_inds <- (nOP+nGP+nEP+1):(nOP+nGP+nEP+nDP)

  parms <- c(occParms,gamParms,epsParms,detParms)
  if(ttdDist == "weibull") parms <- c(parms, "k")
  nP <- length(parms)

  #Likelihood functions--------------------------------------------------------

  nll_R <- function(params){

    #Get occupancy and detection parameters
    psi <- linkFunc(W %*% params[psi_inds])
    psi <- cbind(1-psi, psi)
    lam <- exp(V %*% params[det_inds])

    #Simplified version of Garrard et al. 2013 eqn 5
    #Extended to Weibull
    if(ttdDist=='weibull'){
      k <- exp(params[nP])
      e_lamt <- ( k*lam*(lam*yvec)^(k-1) )^delta * exp(-1*(lam*yvec)^k)
    } else {
      e_lamt <- lam^delta * exp(-lam*yvec)
    }

    #Get probability of y for each z state
    get_Py <- function(e_lamt, delta){
      sum_delt <- as.numeric(sum(delta, na.rm=T)>0)
      c(1-sum_delt, prod(e_lamt[!is.na(e_lamt)]))
    }

    #If dynamic, get col/ext probs and transition prob matrix
    if(T>1){
      col <- plogis(X.gam %*% params[col_inds])
      ext <- plogis(X.eps %*% params[ext_inds])
      phi <- cbind(1-col, col, ext, 1-ext)
    }

    #Begin likelihood calculation
    lik <- rep(NA,N)
    ystart <- 1
    phi_index <- 1
    for (n in 1:N){

      phi_prod <- diag(2)
      #If dynamic, iterate through primary periods
      if(T>1){
        for (t in 1:(T-1)){
          yend <- ystart+J-1
          D_pt <- diag(get_Py(e_lamt[ystart:yend], delta[ystart:yend]))
          phi_t <- matrix(phi[phi_index,],nrow=2, byrow=TRUE)
          phi_prod <- phi_prod %*% ( D_pt %*% phi_t )
          ystart <- ystart + J
          phi_index <- phi_index + 1
        }
      }

      yend <- ystart+J-1
      p_T <- get_Py(e_lamt[ystart:yend], delta[ystart:yend])
      ystart <- ystart + J

      lik[n] <- psi[n,] %*% phi_prod %*% p_T
    }
    -sum(log(lik))
  }

  nll_C <- function(params){
    nll_occuTTD(
          params, yvec, delta, W, V, X.gam, X.eps,
          range(psi_inds)-1, range(det_inds)-1,
          range(col_inds)-1, range(ext_inds)-1,
          linkPsi, ttdDist, N, T, J, naflag
          )
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

  psi <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                          estimates = ests[psi_inds],
                          covMat = as.matrix(covMat[psi_inds,psi_inds]),
                          invlink = invlink,
                          invlinkGrad = linkGrad)

  det <- unmarkedEstimate(name = "Detection", short.name = "lam",
                          estimates = ests[det_inds],
                          covMat = as.matrix(covMat[det_inds,det_inds]),
                          invlink = "exp",
                          invlinkGrad = "exp")

  if(T>1){
    col <- unmarkedEstimate(name = "Colonization", short.name = "col",
                          estimates = ests[col_inds],
                          covMat = as.matrix(covMat[col_inds,col_inds]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")

    ext <- unmarkedEstimate(name = "Extinction", short.name = "ext",
                          estimates = ests[ext_inds],
                          covMat = as.matrix(covMat[ext_inds,ext_inds]),
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad")


    estimateList <- unmarkedEstimateList(list(psi = psi, col = col,
                                            ext = ext, det=det))
  } else {
    estimateList <- unmarkedEstimateList(list(psi = psi, det=det))
  }

  #Add Weibull shape parameter if necessary
  if(ttdDist=="weibull"){
    estimateList@estimates$shape <- unmarkedEstimate(name = "Weibull shape",
    short.name = "k", estimates = ests[nP],
    covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
    invlinkGrad = "exp")
  }

  umfit <- new("unmarkedFitOccuTTD", fitType = "occuTTD",
               call = match.call(),
               formula = formula,
               psiformula = psiformula,
               gamformula = gammaformula,
               epsformula = epsilonformula,
               detformula = detformula,
               data = data, sitesRemoved = removed,
               estimates = estimateList,
               AIC = fmAIC, opt = fm, negLogLike = fm$value,
               nllFun = nll)

  return(umfit)
}
