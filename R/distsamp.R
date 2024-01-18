
distsamp <- function(formula, data,
    keyfun=c("halfnorm", "exp", "hazard", "uniform"),
    output=c("density", "abund"), unitsOut=c("ha", "kmsq"), starts,
    method="BFGS", se = TRUE, engine = c("C", "R", "TMB"),
    rel.tol=0.001, ...)
{

    # Check arguments
    engine <- match.arg(engine)
    if(any(sapply(split_formula(formula), has_random))) engine <- "TMB"
    keyfun <- match.arg(keyfun)
    output <- match.arg(output)
    unitsOut <- match.arg(unitsOut)

    if(missing(starts)) starts <- NULL

    #Generate design matrix
    designMats <- getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if(is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    if(is.null(V.offset))
        V.offset <- rep(0, nrow(V))

    M <- nrow(y)
    J <- ncol(y)

    # Distance sampling design info
    db <- data@dist.breaks
    tlength <- data@tlength
    survey <- data@survey
    w <- diff(db)
    unitsIn <- data@unitsIn
    u <- a <- matrix(NA, M, J)
    switch(survey,
        line = {
            for(i in 1:M) {
                a[i,] <- tlength[i] * w
                u[i,] <- a[i,] / sum(a[i,])
                }
            },
        point = {
            for(i in 1:M) {
                a[i, 1] <- pi*db[2]^2
                for(j in 2:J)
                    a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
                u[i,] <- a[i,] / sum(a[i,])
                }
            })
    switch(survey,
        line = A <- rowSums(a) * 2,
        point = A <- rowSums(a))
    switch(unitsIn,
        m = A <- A / 1e6,
        km = A <- A)
    switch(unitsOut,
        ha = A <- A * 100,
        kmsq = A <- A)

    # Set up parameters
    lamParms <- colnames(X)
    detParms <- colnames(V)
    scaleParms <- character(0)
    nAP <- length(lamParms)
    nDP <- length(detParms)
    nP <- nAP + nDP
    lamIdx <- 1:nAP
    detIdx <- (nAP+1):nP
    starts_default <- c(rep(0, nAP), log(max(db)), rep(0, nDP-1))

    if(keyfun=="uniform"){
      detParms <- character(0)
      detIdx <- numeric(0)
      starts_default <- rep(0, nAP)
    }
    if(keyfun=="exp"){
      starts_default[(nAP+1)] <- 0
      # maybe this should be default everywhere
      if(engine == "TMB") starts_default[(nAP+1)] <- log(median(db))
    }
    if(keyfun=="hazard"){
      nP <- nP + 1
      scaleParms <- "scale"
      starts_default[(nAP+1)] <- log(median(db))
      starts_default <- c(starts_default, 1)
    }
    names(starts_default) <- c(lamParms, detParms, scaleParms)

    if(engine=="R") {

    cp <- matrix(NA, M, J)
    switch(keyfun,
    halfnorm = {
        nll <- function(param) {
            sigma <- drop(exp(V %*% param[(nAP+1):nP] + V.offset))
            lambda <- drop(exp(X %*% param[1:nAP] + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    f.0 <- 2 * dnorm(0, 0, sd=sigma[i])
                    int <- 2 * (pnorm(db[-1], 0, sd=sigma[i]) -
                        pnorm(db[-(J+1)], 0, sd=sigma[i]))
                    cp[i,] <- int / f.0 / w
                    },
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhn, db[j], db[j+1], sigma=sigma[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    exp = {
        nll <- function(param) {
            rate <- drop(exp(V %*% param[(nAP+1):nP] + V.offset))
            lambda <- drop(exp(X %*% param[1:nAP] + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxexp, db[j], db[j+1], rate=rate[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value / w[j]
                        else {
                            cp[i, j] <- NA
                            }
                        }},
                point = {
                    for(j in 1:J) {
                        int <- integrate(grexp, db[j], db[j+1], rate=rate[i],
                            stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    hazard = {
        nll <- function(param) {
            shape <- drop(exp(V %*% param[(nAP+1):(nP-1)] + V.offset))
            scale <- drop(exp(param[nP]))
            lambda <- drop(exp(X %*% param[1:nAP] + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        int <- integrate(gxhaz, db[j], db[j+1], shape=shape[i],
                            scale=scale, stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value / w[j]
                        else {
                            cp[i, j] <- NA
                            }
                        }},
                point = {
                    for(j in 1:J) {
                        int <- integrate(grhaz, db[j], db[j+1], shape=shape[i],
                            scale=scale, stop.on.error=FALSE)
                        mess <- int$message
                        if(identical(mess, "OK"))
                            cp[i, j] <- int$value * 2*pi / a[i,j]
                        else {
                            cp[i, j] <- NA
                            }
                        }

                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            ll <- dpois(y, lambda * cp, log=TRUE)
            -sum(ll)
            }},
    uniform = {
        nll <- function(param) {
            lambda <- drop(exp(X %*% param + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            ll <- dpois(y, lambda * u, log=TRUE)
            -sum(ll)
            }
        })
    } else if(engine=="C") {
        nll <- function(param) {
            beta.lam <- param[1:nAP]
            if(identical(keyfun, "hazard")) {
                beta.sig <- param[(nAP+1):(nP-1)]
                scale <- exp(param[nP])
            } else {
                beta.sig <- param[(nAP+1):nP]
                scale <- -99.0
            }
            lambda <- drop(exp(X %*% beta.lam + X.offset))
            if(identical(output, "density"))
                lambda <- lambda * A
            sigma <- drop(exp(V %*% beta.sig + V.offset))
            nll_distsamp(
                  y, lambda, sigma, scale,
                  a, u, w, db,
                  keyfun, survey
            )
        }
    }

    if(engine %in% c("C","R")){
      if(is.null(starts)) starts <- starts_default
      fm <- optim(starts, nll, method=method, hessian=se, ...)

      ests <- fm$par
      names(ests) <- c(lamParms, detParms, scaleParms)
      covMat <- invertHessian(fm, nP, se)
      fmAIC <- 2 * fm$value + 2 * nP
      tmb_mod <- NULL

      # Organize fixed-effect estimates
      state_coef <- list(ests=ests[lamIdx], cov=as.matrix(covMat[lamIdx, lamIdx]))

      if(keyfun != "uniform"){
        det_coef <- list(ests=ests[detIdx], cov=as.matrix(covMat[detIdx, detIdx]))
      }
      if(keyfun == "hazard") {
        scale_coef <- list(ests=ests[nP], cov=as.matrix(covMat[nP,nP]))
      }

      # No random effects in C or R engines
      state_rand_info <- det_rand_info <- list()

    } else if(engine == "TMB"){

      # Set up TMB input data
      if(output == "abund") A <- rep(1, length(A))
      forms <- split_formula(formula)
      inps <- get_ranef_inputs(forms, list(det=siteCovs(data), state=siteCovs(data)),
                               list(V, X), designMats[c("Z_det","Z_state")])

      keyfun_type <- switch(keyfun, uniform={0}, halfnorm={1}, exp={2},
                            hazard={3})
      survey_type <- switch(survey, line={0}, point={1})
      tmb_dat <- c(list(y=y, survey_type=survey_type, keyfun_type=keyfun_type,
                        A=A, db=db, a=a, w=w, u=u, offset_state=X.offset,
                        offset_det=V.offset), inps$data)

      tmb_param <- c(inps$pars, list(beta_scale=rep(0,0)))

      if(is.null(starts)){
        if(keyfun != "uniform") tmb_param$beta_det[1] <- log(median(db))
      }

      if(keyfun == "hazard") tmb_param$beta_scale <- rep(0,1)
      if(keyfun == "uniform") tmb_param$beta_det <- rep(0,0)

      # Fit model in TMB
      tmb_out <- fit_TMB("tmb_distsamp", tmb_dat, tmb_param, inps$rand_ef,
                         starts=starts, method)
      tmb_mod <- tmb_out$TMB
      fm <- tmb_out$opt
      fmAIC <- tmb_out$AIC
      nll <- tmb_mod$fn

      # Organize fixed-effect estimate from TMB output
      state_coef <- get_coef_info(tmb_out$sdr, "state", lamParms, lamIdx)
      det_coef <- get_coef_info(tmb_out$sdr, "det", detParms, detIdx)

      if(keyfun=="hazard"){
        scale_coef <- get_coef_info(tmb_out$sdr, "scale", scaleParms, nP)
      }

      # Organize random-effect estimates from TMB output
      state_rand_info <- get_randvar_info(tmb_out$sdr, "state", forms[[2]], siteCovs(data))
      det_rand_info <- get_randvar_info(tmb_out$sdr, "det", forms[[1]], siteCovs(data))

    }

    stateName <- switch(output, abund = "Abundance", density = "Density")
    stateEstimates <- unmarkedEstimate(name = stateName,
        short.name = "lam", estimates = state_coef$ests, covMat = state_coef$cov,
        fixed=1:nAP, invlink = "exp", invlinkGrad = "exp", randomVarInfo=state_rand_info)
    estimateList <- unmarkedEstimateList(list(state=stateEstimates))

    if(keyfun != "uniform") {
      detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
        estimates = det_coef$ests, covMat = det_coef$cov, fixed=1:nDP,
        invlink = "exp", invlinkGrad = "exp", randomVarInfo=det_rand_info)
      estimateList@estimates$det <- detEstimates

      if(keyfun == "hazard"){
        scaleEstimates <- unmarkedEstimate(name = "Hazard-rate(scale)",
          short.name = "p", estimates = scale_coef$ests,
          covMat = scale_coef$cov, fixed=1, invlink = "exp", invlinkGrad = "exp",
          randomVarInfo=list())
        estimateList@estimates$scale <- scaleEstimates
      }
    }

    dsfit <- new("unmarkedFitDS", fitType = "distsamp", call = match.call(),
        opt = fm, formula = formula, data = data, keyfun=keyfun,
        sitesRemoved = designMats$removed.sites, unitsOut=unitsOut,
        estimates = estimateList, AIC = fmAIC, negLogLike = fm$value,
        nllFun = nll, output=output, TMB=tmb_mod)
    return(dsfit)
}


# Detection functions

gxhn <- function(x, sigma) exp(-x^2/(2 * sigma^2))
gxexp <- function(x, rate) exp(-x / rate)
gxhaz <- function(x, shape, scale)  1 - exp(-(x/shape)^-scale)
grhn <- function(r, sigma) exp(-r^2/(2 * sigma^2)) * r
grexp <- function(r, rate) exp(-r / rate) * r
grhaz <- function(r, shape, scale)  (1 - exp(-(r/shape)^-scale)) * r

dxhn <- function(x, sigma)
	gxhn(x=x, sigma=sigma) / integrate(gxhn, 0, Inf, sigma=sigma)$value
drhn <- function(r, sigma)
	grhn(r=r, sigma=sigma) / integrate(grhn, 0, Inf, sigma=sigma)$value
dxexp <- function(x, rate)
	gxexp(x=x, rate=rate) / integrate(gxexp, 0, Inf, rate=rate)$value
drexp <- function(r, rate)
	grexp(r=r, rate=rate) / integrate(grexp, 0, Inf, rate=rate)$value
dxhaz <- function(x, shape, scale)
	gxhaz(x=x, shape=shape, scale=scale) / integrate(gxhaz, 0, Inf,
		shape=shape, scale=scale)$value
drhaz <- function(r, shape, scale)
	grhaz(r=r, shape=shape, scale=scale) / integrate(grhaz, 0, Inf,
		shape=shape, scale=scale)$value


