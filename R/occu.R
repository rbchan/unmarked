
#  Fit the occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0),
                 linkPsi = c("logit", "cloglog"), starts, method = "BFGS",
                 se = TRUE, engine = c("C", "R", "TMB", "RTMB"), threads=1, ...) {

  # Check arguments------------------------------------------------------------
  if(!is(data, "unmarkedFrameOccu"))
    stop("Data is not an unmarkedFrameOccu object.")

  engine <- match.arg(engine, c("C", "R", "TMB", "RTMB"))
  if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
  if(length(knownOcc)>0 & engine == "TMB"){
    stop("TMB engine does not support knownOcc argument", call.=FALSE)
  }

  linkPsi <- match.arg(linkPsi, c("logit","cloglog"))
  psiLinkFunc <- ifelse(linkPsi=="cloglog", cloglog, plogis)
  psiInvLink <- ifelse(linkPsi=="cloglog", "cloglog", "logistic")
  psiLinkGrad <- ifelse(linkPsi=="cloglog", "cloglog.grad", "logistic.grad")

  # Format input data----------------------------------------------------------
  designMats <- getDesign(data, formula)
  y <- truncateToBinary(designMats$y)
  X <- designMats$X; V <- designMats$V;
  Z_state <- designMats$Z_state; Z_det <- designMats$Z_det
  removed <- designMats$removed.sites
  X.offset <- designMats$X.offset; V.offset <- designMats$V.offset

  # Re-format some variables for C and R engines
  yvec <- as.numeric(t(y))
  navec <- is.na(yvec)
  nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i

  # convert knownOcc to logical so we can correctly to handle NAs.
  knownOccLog <- rep(FALSE, numSites(data))
  knownOccLog[knownOcc] <- TRUE
  if(length(removed)>0) knownOccLog <- knownOccLog[-removed]

  # Set up parameter names and indices-----------------------------------------
  occParms <- colnames(X)
  detParms <- colnames(V)
  nDP <- ncol(V)
  nOP <- ncol(X)
  nP <- nDP + nOP
  psiIdx <- 1:nOP
  pIdx <- (nOP+1):nP

  # Set up negative log likelihood functions for C++ and R engines-----_-------
  if(identical(engine, "C")) {
    nll <- function(params) {
      beta.psi <- params[1:nOP]
      beta.p <- params[(nOP+1):nP]
      nll_occu(
        yvec, X, V, beta.psi, beta.p, nd, knownOccLog, navec,
        X.offset, V.offset, linkPsi
      )
    }
  } else if (identical(engine, "R")){

    J <- ncol(y)
    M <- nrow(y)

    nll <- function(params) {
      psi <- psiLinkFunc(X %*% params[1 : nOP] + X.offset)
      psi[knownOccLog] <- 1
      pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
      cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
      cp[navec] <- 1 # so that NA's don't modify likelihood
      cpmat <- matrix(cp, M, J, byrow = TRUE) #
      loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
      -sum(loglik)
    }
  }

  # Fit model with C++ and R engines-------------------------------------------
  if(engine %in% c("C", "R")){
    if(missing(starts)) starts <- rep(0, nP)
    if(length(starts) != nP){
      stop(paste("The number of starting values should be", nP))
    }
    fm <- optim(starts, nll, method = method, hessian = se, ...)
    covMat <- invertHessian(fm, nP, se)
    ests <- fm$par
    tmb_mod <- NULL
    fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)

    # Organize effect estimates
    names(ests) <- c(occParms, detParms)
    state_coef <- list(ests=ests[1:nOP], cov=as.matrix(covMat[1:nOP,1:nOP]))
    det_coef <- list(ests=ests[(nOP+1):nP],
                     cov=as.matrix(covMat[(nOP+1):nP, (nOP+1):nP]))

    # No random effects for C++ and R engines
    state_rand_info <- det_rand_info <- list()

  # Fit model with TMB engine--------------------------------------------------
  } else if(identical(engine, "TMB")){

    # Set up TMB input data
    forms <- split_formula(formula)
    obs_all <- add_covariates(obsCovs(data), siteCovs(data), length(getY(data)))
    inps <- get_ranef_inputs(forms, list(det=obs_all, state=siteCovs(data)),
                             list(V, X), designMats[c("Z_det","Z_state")])

    tmb_dat <- c(list(y=y, no_detect=nd, link=ifelse(linkPsi=="cloglog",1,0),
                      offset_state=X.offset, offset_det=V.offset), inps$data)

    # Fit model with TMB
    if(missing(starts)) starts <- NULL
    tmb_out <- fit_TMB("tmb_occu", tmb_dat, inps$pars, inps$rand_ef,
                        starts=starts, method, ...)
    tmb_mod <- tmb_out$TMB
    fm <- tmb_out$opt
    fmAIC <- tmb_out$AIC
    nll <- tmb_mod$fn

    # Organize fixed effect estimates
    state_coef <- get_coef_info(tmb_out$sdr, "state", occParms, 1:nOP)
    det_coef <- get_coef_info(tmb_out$sdr, "det", detParms, (nOP+1):nP)

    # Organize random effect estimates
    state_rand_info <- get_randvar_info(tmb_out$sdr, "state", forms[[2]],
                                        siteCovs(data))
    det_rand_info <- get_randvar_info(tmb_out$sdr, "det", forms[[1]],
                                      obs_all)

  } else if(identical(engine, "RTMB")){

    # Set up TMB input data
    forms <- split_formula(formula)
    obs_all <- add_covariates(obsCovs(data), siteCovs(data), length(getY(data)))
    inps <- get_ranef_inputs(forms, list(det=obs_all, state=siteCovs(data)),
                             list(V, X), designMats[c("Z_det","Z_state")])

    tmb_dat <- c(list(y=y, no_detect=nd, link=ifelse(linkPsi=="cloglog",1,0),
                      offset_state=X.offset, offset_det=V.offset), inps$data)

    f <- function(pars){
      beta_state <- pars$beta_state
      beta_det <- pars$beta_det

      M <- nrow(tmb_dat$y)
      J <- ncol(tmb_dat$y)

      psi <- tmb_dat$X_state %*% beta_state + tmb_dat$offset_state
      psi <- RTMB::plogis(psi)
      #psi <- 1/(1+exp(-psi))  

      p <- tmb_dat$X_det %*% beta_det + tmb_dat$offset_det
      #p <- 1/(1+exp(-p))
      p <- RTMB::plogis(p)

      nll <- 0

      for (i in 1:M){
        cp <- 1.0
        pind <- (i-1) * J + 1
        #for (j in 1:J){
        #  cp <- cp * RTMB::dbinom(tmb_dat$y[i,j], 1, p[pind])
        #  pind <- pind + 1
        #}
        cp <- prod(RTMB::dbinom(tmb_dat$y[i,], 1, p[pind:(pind+J-1)]))
        nll <- nll - log(psi[i] * cp + (1-psi[i]) * tmb_dat$no_detect[i])
      }

      nll
    }

    tmb_mod <- RTMB::MakeADFun(f, parameters=inps$pars, silent=TRUE)
    opt <- optim(tmb_mod$par, fn=tmb_mod$fn, gr=tmb_mod$gr, method=method, ...)
    sdr <- RTMB::sdreport(tmb_mod, getJointPrecision=TRUE)
    sdr$par <- tmb_mod$par

    fm <- opt
    fmAIC <- 2 * opt$value + 2 * length(unlist(inps$pars)) # wrong?
    nll <- tmb_mod$fn

    # Organize fixed effect estimates
    state_coef <- get_coef_info(sdr, "state", occParms, 1:nOP)
    det_coef <- get_coef_info(sdr, "det", detParms, (nOP+1):nP)

    # Organize random effect estimates
    state_rand_info <- get_randvar_info(sdr, "state", forms[[2]],
                                        siteCovs(data))
    det_rand_info <- get_randvar_info(sdr, "det", forms[[1]],
                                      obs_all)

  }

  # Create unmarkedEstimates---------------------------------------------------
  state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
                            estimates = state_coef$est,
                            covMat = state_coef$cov,
                            fixed = 1:nOP,
                            invlink = psiInvLink,
                            invlinkGrad = psiLinkGrad,
                            randomVarInfo=state_rand_info
                            )

  det <- unmarkedEstimate(name = "Detection", short.name = "p",
                          estimates =det_coef$est,
                          covMat = det_coef$cov,
                          fixed = 1:nDP,
                          invlink = "logistic",
                          invlinkGrad = "logistic.grad",
                          randomVarInfo=det_rand_info
                          )

  estimateList <- unmarkedEstimateList(list(state=state, det=det))

  # Create unmarkedFit object--------------------------------------------------
  umfit <- new("unmarkedFitOccu", fitType = "occu", call = match.call(),
                 formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = fm,
                 negLogLike = fm$value,
                 nllFun = nll, knownOcc = knownOccLog, TMB=tmb_mod)

  return(umfit)
}
