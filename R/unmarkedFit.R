setClass("unmarkedFit",
   representation(fitType = "character",
        call = "call",
        formula = "formula",
        data = "unmarkedFrame",
        sitesRemoved = "numeric",  # vector of indices of removed sites
        estimates = "unmarkedEstimateList",
        AIC = "numeric",
        opt = "list",
        negLogLike = "numeric",
        nllFun = "function",
        bootstrapSamples = "optionalList",
        covMatBS = "optionalMatrix", # list of bootstrap sample fits
        TMB = "optionalList")) #TMB output object

# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call, formula, data, sitesRemoved,
    estimates, AIC, opt, negLogLike, nllFun)
{
    umfit <- new("unmarkedFit", fitType = fitType, call = call,
        formula = formula, data = data, sitesRemoved = sitesRemoved,
        estimates = estimates, AIC = AIC, opt = opt,
        negLogLike = negLogLike,
        nllFun = nllFun)
    return(umfit)
}

# ---------------------------- CHILD CLASSES ----------------------------


setClass("unmarkedFitDS",
    representation(
        keyfun = "character",
        unitsOut = "character",
        output = "character"),
        contains = "unmarkedFit")



setClass("unmarkedFitPCount",
    representation(
        K = "numeric",
        mixture = "character"),
    contains = "unmarkedFit")

# This class is not used directly, just used as a base for for PCO, MMO, DSO
setClass("unmarkedFitDailMadsen",
        representation(
            K = "numeric",
            mixture = "character",
            formlist = "list",
            dynamics = "character",
            immigration = "logical",
            fix = "character"),
         contains = "unmarkedFit")

setClass("unmarkedFitPCO", contains = "unmarkedFitDailMadsen")

setClass("unmarkedFitMMO", contains = "unmarkedFitDailMadsen")

setClass("unmarkedFitDSO",
        representation(
            keyfun = "character",
            unitsOut = "character",
            output = "character"),
        contains = "unmarkedFitDailMadsen")

setClass("unmarkedFitOccu",
    representation(knownOcc = "logical"),
    contains = "unmarkedFit")

setClass("unmarkedFitOccuPEN",
    representation(
	knownOcc = "logical",
	pen.type = "character",
	lambda = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitOccuPEN_CV",
    representation(
	knownOcc = "logical",
	pen.type = "character",
	lambdaVec = "numeric",
	k = "numeric",
	foldAssignments = "numeric",
	lambdaScores = "numeric",
	chosenLambda = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitOccuFP",
         representation(knownOcc = "logical",
            detformula = "formula",
            FPformula = "formula",
            Bformula = "formula",
            stateformula = "formula",
            type = "numeric"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuMulti",
         representation(
            detformulas = "character",
            stateformulas = "character"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuMS",
         representation(
            detformulas = "character",
            psiformulas = "character",
            phiformulas = "character",
            parameterization = "character"),
         contains = "unmarkedFit")

setClass("unmarkedFitOccuTTD",
    representation(
        psiformula = "formula",
        gamformula = "formula",
        epsformula = "formula",
        detformula = "formula"),
    contains = "unmarkedFit")

setClass("unmarkedFitNmixTTD",
         representation(
           stateformula = "formula",
           detformula = "formula",
           K = "numeric"),
         contains = "unmarkedFit")

setClass("unmarkedFitMPois",
    contains = "unmarkedFit")


setClass("unmarkedFitOccuRN",
    representation(
      K = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitMNmix",
    representation(constraint = "numeric"),
    contains = "unmarkedFit")

setClass("unmarkedFitColExt",
    representation(
        phi = "matrix",
        psiformula = "formula",
        gamformula = "formula",
        epsformula = "formula",
        detformula = "formula",
        projected = "array",
        projected.mean = "matrix",
        smoothed = "array",
        smoothed.mean = "matrix",
        projected.mean.bsse = "optionalMatrix",
        smoothed.mean.bsse = "optionalMatrix"),
    contains = "unmarkedFit")


setClass("unmarkedFitGMM",
    representation(
        formlist = "list",
        mixture = "character",
        K = "numeric"),
    contains = "unmarkedFit")


setClass("unmarkedFitGDS",
    representation(
        keyfun = "character",
        unitsOut = "character",
        output = "character"),
    contains = "unmarkedFitGMM")


setClass("unmarkedFitGPC",
    contains = "unmarkedFitGMM")



# -------------------------- Show and Summary ----------------------------


setMethod("show", "unmarkedFit", function(object)
{
    cat("\nCall:\n")
    print(object@call)
    cat("\n")
    show(object@estimates)
    cat("AIC:", object@AIC,"\n")
    if(!identical(object@opt$convergence, 0L))
        warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
})




setMethod("summary", "unmarkedFit", function(object)
{
    cat("\nCall:\n")
    print(object@call)
    cat("\n")
    summaryOut <- summary(object@estimates)
    cat("AIC:", object@AIC,"\n")
    cat("Number of sites:", sampleSize(object))
    if(length(object@sitesRemoved) > 0)
        cat("\nID of sites removed due to NA:", object@sitesRemoved)
    cat("\noptim convergence code:", object@opt$convergence)
    cat("\noptim iterations:", object@opt$counts[1], "\n")
    if(!identical(object@opt$convergence, 0L))
    warning("Model did not converge. Try providing starting values or increasing maxit control argment.")
    cat("Bootstrap iterations:", length(object@bootstrapSamples), "\n\n")
    invisible(summaryOut)
})



setMethod("summary", "unmarkedFitDS", function(object)
{
    callNextMethod()
    cat("Survey design: ", object@data@survey, "-transect", sep="")
    cat("\nDetection function:", object@keyfun)
    cat("\nUnitsIn:", object@data@unitsIn)
    cat("\nUnitsOut:", object@unitsOut, "\n\n")
})




# Compute linear combinations of estimates in unmarkedFit objects.
setMethod("linearComb",
    signature(obj = "unmarkedFit", coefficients = "matrixOrVector"),
    function(obj, coefficients, type, offset = NULL, re.form=NULL)
{
    stopifnot(!missing(type))
    stopifnot(type %in% names(obj))
    estimate <- obj@estimates[type]
    linearComb(estimate, coefficients, offset, re.form)
})


setMethod("backTransform", "unmarkedFit", function(obj, type)
{
    est <- obj[type]
    if(length(est@estimates) == 1) {
        lc <- linearComb(est, 1)
        return(backTransform(lc))
    } else {
        stop('Cannot directly backTransform an unmarkedEstimate with length > 1.')
        }
})


setMethod("[", "unmarkedFit",
          function(x, i, j, drop) {
              x@estimates[i]
          })


setMethod("names", "unmarkedFit",
          function(x) {
              names(x@estimates)
          })



# ---------------------- coef, vcov, and SE ------------------------------


setMethod("coef", "unmarkedFit",
    function(object, type, altNames = TRUE, fixedOnly=TRUE)
{
    if(missing(type)) {
        co <- lapply(object@estimates@estimates,
            function(x) coef(x, altNames=altNames, fixedOnly=fixedOnly))
        names(co) <- NULL
        co <- unlist(co)
    } else {
        co <- coef(object[type], altNames=altNames, fixedOnly=fixedOnly)
        }
    co
})


setMethod("vcov", "unmarkedFit",
    function (object, type, altNames = TRUE, method = "hessian", fixedOnly=TRUE, ...)
{
    method <- match.arg(method, c("hessian", "nonparboot"))
    if(.hasSlot(object, "TMB") && !is.null(object@TMB)) method <- "TMB"
    switch(method,
           hessian = {
            if (is.null(object@opt$hessian)) {
                stop("Hessian was not computed for this model.")
            }
            v <- solve(hessian(object))
        },
        nonparboot = {
            if (is.null(object@bootstrapSamples)) {
                stop("No bootstrap samples have been drawn. Use nonparboot first.")
                }
            v <- object@covMatBS
        },
        TMB = {
          return(vcov_TMB(object, type, fixedOnly))
        })
    rownames(v) <- colnames(v) <- names(coef(object, altNames=altNames))
    if (missing(type)) {
        return (v)
    } else {
        inds <- .estimateInds(object)[[type]]
        return (v[inds, inds, drop = FALSE])
        }
})

setMethod("vcov", "unmarkedFitOccuMulti",
    function (object, type, altNames = TRUE, method = "hessian", ...)
{

  pen <- object@call[["penalty"]]
  if(is.null(pen)) pen <- 0
  if(pen>0) method <- "nonparboot"
  callNextMethod(object, type, altNames, method=method, ...)

})


setMethod("SE", "unmarkedFit", function(obj,...)
{
    v <- vcov(obj,...)
    sqrt(diag(v))
})



setMethod("logLik", "unmarkedFit", function(object, ...)
{
    if(length(list(...)))
        warning("extra arguments discarded")
    ll <- -object@negLogLike
    #attr(ll, "df") <- length(coef(object))
    #class(ll) <- "logLik"
    return(ll)
})



setMethod("LRT", c(m1="unmarkedFit", m2="unmarkedFit"), function(m1, m2)
{
    ll1 <- logLik(m1)
    ll2 <- logLik(m2)
    chisq <- 2 * abs(ll1 - ll2)
    DF <- abs(length(coef(m1)) - length(coef(m2)))
    pval <- pchisq(chisq, DF, lower.tail=FALSE)
    return(data.frame(Chisq=chisq, DF = DF, 'Pr(>Chisq)' = pval,
        check.names=F))
})





setMethod("confint", "unmarkedFit", function(object, parm, level = 0.95,
    type, method = c("normal", "profile"))
{
    method <- match.arg(method)
    if(missing(type))
        stop(paste("Must specify type as one of (", paste(names(object), collapse=", "),").",sep=""))
    if(missing(parm))
        parm <- object[type]@fixed
    if(method == "normal") {
        callGeneric(object[type],parm = parm, level = level)
    } else {
        nllFun <- nllFun(object)
        ests <- mle(object)
        nP <- length(parm)
        ci <- matrix(NA, nP, 2)

        ## create table to match parm numbers with est/parm numbers
        types <- names(object)
        numbertable <- data.frame(type = character(0), num = numeric(0))
        for(i in seq(length=length(types))) {
            length.est <- length(object[i]@estimates)
            numbertable <- rbind(numbertable, data.frame(type =
                rep(types[i], length.est), num = seq(length=length.est)))
            }
        parm.fullnums <- which(numbertable$type == type &
            numbertable$num %in% parm)

        for(i in seq(length=nP)) {
            cat("Profiling parameter",i,"of",nP,"...")
            se <- SE(object[type])
            whichPar <- parm.fullnums[i]
            ci[i,] <- profileCI(nllFun, whichPar=whichPar, MLE=ests,
                interval=ests[whichPar] + 10*se[i]*c(-1,1), level=level)
            cat(" done.\n")
            }
        rownames(ci) <- names(coef(object[type]))[parm]
        colnames(ci) <- c((1-level)/2, 1- (1-level)/2)
        return(ci)
        }
})



setMethod("fitted", "unmarkedFit",
    function(object, na.rm = FALSE)
{
    data <- object@data
    des <- getDesign(data, object@formula, na.rm = na.rm)
    X <- cbind(des$X, des$Z_state)
    X.offset <- des$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    beta_state <- coef(object, 'state', fixedOnly=FALSE)
    state <- do.call(object['state']@invlink,
        list(as.matrix(X %*% beta_state + X.offset)))
    state <- as.numeric(state)  ## E(X) for most models
    p <- getP(object, na.rm = na.rm) # P(detection | presence)
    fitted <- state * p  # true for models with E[Y] = p * E[X]
    fitted
})



setMethod("fitted", "unmarkedFitOccuFP", function(object, na.rm = FALSE)
{
  cat("fitted is not implemented for occuFP at this time")
})


setMethod("fitted", "unmarkedFitDS", function(object, na.rm = FALSE)
{
    data <- object@data
    db <- data@dist.breaks
    w <- diff(db)
    D <- getDesign(data, object@formula, na.rm = na.rm)
    X <- D$X
    X.offset <- D$X.offset
    if (is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    M <- nrow(X)
    J <- length(w)
    lambda <- drop(exp(X %*% coef(object, 'state') + X.offset))
    if(identical(object@output, "density")) {
        a <- matrix(NA, M, J)
        switch(data@survey,
            line = {
                tlength <- data@tlength
                A <- tlength * max(db) * 2
                },
            point = {
                A <- pi * max(db)^2
                })
        switch(data@unitsIn,
            m = A <- A / 1e6,
            km = A <- A)
        switch(object@unitsOut,
            m = A <- A * 1e6,
            ha = A <- A * 100,
            kmsq = A <- A)
        lambda <- lambda * A
        }
    cp <- getP(object, na.rm = na.rm)
    fitted <- lambda * cp
    fitted
})



setMethod("fitted", "unmarkedFitOccu", function(object, na.rm = FALSE)
{
    data <- object@data
    des <- getDesign(data, object@formula, na.rm = na.rm)
    X <- des$X
    X.offset <- des$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    state <- plogis(X %*% coef(object, 'state') + X.offset)
    state <- as.numeric(state)  ## E(X) for most models
    state[object@knownOcc] <- 1
    p <- getP(object, na.rm = na.rm) # P(detection | presence)
    fitted <- state * p  # true for models with E[Y] = p * E[X]
    fitted
})



setMethod("fitted", "unmarkedFitPCount", function(object, K, na.rm = FALSE)
{
    data <- object@data
    des <- getDesign(data, object@formula, na.rm = na.rm)
    X <- cbind(des$X, des$Z_state)
    X.offset <- des$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    y <- des$y	# getY(data) ... to be consistent w/NA handling?
    M <- nrow(X)
    J <- ncol(y)
    state <- exp(X %*% coef(object, 'state', fixedOnly=FALSE) + X.offset)
    p <- getP(object, na.rm = na.rm)
    mix <- object@mixture
##    if(!is.missing(K))
##        warning("The K argument is ignored")
    switch(mix,
           P = {
               fitted <- as.numeric(state) * p
           },
           NB = {
               ## I don't think this sum is necessary. Could do:
               ## fitted <- as.numeric(state) * p
               if(missing(K))
                   K <- object@K
               k <- 0:K
               k.ijk <- rep(k, M*J)
               state.ijk <- state[rep(1:M, each = J*(K+1))]
               alpha <- exp(coef(object['alpha']))
               prob.ijk <- dnbinom(k.ijk, mu = state.ijk, size = alpha)
               all <- cbind(rep(as.vector(t(p)), each = K + 1), k.ijk,
                            prob.ijk)
               prod.ijk <- rowProds(all)
               fitted <- colSums(matrix(prod.ijk, K + 1, M*J))
               fitted <- matrix(fitted, M, J, byrow = TRUE)
           },
           ZIP = {
               psi <- plogis(coef(object['psi']))
               lambda <- as.numeric(state)
               E.N <- (1-psi)*lambda
#               fitted <- (1-psi)*lambda
#               fitted <- matrix(fitted, M, J, byrow=TRUE) # BUG
               fitted <- E.N * p
           })
    return(fitted)
})


setMethod("fitted", "unmarkedFitDailMadsen",
    function(object, K, na.rm = FALSE)
{
    dynamics <- object@dynamics
    mixture <- object@mixture
    #To partially handle old saved model objects
    fix <- tryCatch(object@fix, error=function(e) "none")
    immigration <- tryCatch(object@immigration, error=function(e) FALSE)
    data <- getData(object)
    D <- getDesign(data, object@formula, na.rm = na.rm)
    Xlam <- D$Xlam; Xgam <- D$Xgam; Xom <- D$Xom; Xiota <- D$Xiota
    Xlam.offset <- D$Xlam.offset; Xgam.offset <- D$Xgam.offset
    Xom.offset <- D$Xom.offset
    Xiota.offset <- D$Xiota.offset
    delta <- D$delta #FIXME this isn't returned propertly when na.rm=F

    y <- D$y
    M <- nrow(y)
    T <- data@numPrimary
    J <- ncol(y) / T

    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
    if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
    if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
    if(is.null(Xiota.offset)) Xiota.offset <- rep(0, M*(T-1))

    lambda <- exp(Xlam %*% coef(object, 'lambda') + Xlam.offset)
    if(identical(mixture, "ZIP")) {
        psi <- plogis(coef(object, type="psi"))
        lambda <- (1-psi)*lambda
    }
    if (fix == 'omega'){
      omega <- matrix(1, M, T-1)
    } else if(!identical(dynamics, "trend")) {
        if(identical(dynamics, "ricker") || identical(dynamics, "gompertz"))
            omega <- matrix(exp(Xom %*% coef(object, 'omega') + Xom.offset),
                        M, T-1, byrow=TRUE)
        else
            omega <- matrix(plogis(Xom %*% coef(object, 'omega') + Xom.offset),
                        M, T-1, byrow=TRUE)
    }
    if(fix == "gamma"){
        gamma <- matrix(0, M, T-1)
    } else if(!identical(dynamics, "notrend")){
        gamma <- matrix(exp(Xgam %*% coef(object, 'gamma') + Xgam.offset),
                        M, T-1, byrow=TRUE)
    } else {
        if(identical(dynamics, "notrend"))
            gamma <- (1-omega)*lambda
        }
    if(immigration)
        iota <- matrix(exp(Xiota %*% coef(object, 'iota') + Xiota.offset),
                        M, T-1, byrow=TRUE)
    else
        iota <- matrix(0, M, T-1)

    N <- matrix(NA, M, T)
    for(i in 1:M) {
        N[i, 1] <- lambda[i]
        if(delta[i, 1] > 1) {
            for(d in 2:delta[i ,1]) {
                if(identical(dynamics, "autoreg"))
                    N[i, 1] <- N[i, 1] * (omega[i,1] + gamma[i, 1]) + iota[i, 1]
            else if(identical(dynamics, "trend"))
                N[i,1] <- N[i,1] * gamma[i,1] + iota[i, 1]
            else if(identical(dynamics, "ricker"))
                N[i,1] <- N[i,1] * exp(gamma[i,1]*(1-N[i,1]/omega[i,1])) +
                    iota[i, 1]
            else if(identical(dynamics, "gompertz"))
                N[i,1] <- N[i,1] * exp(gamma[i,1]*(1-log(N[i,1]+1)/
                  log(omega[i,1]+1))) + iota[i, 1]
            else
                N[i,1] <- N[i,1] * omega[i,1] + gamma[i,1]
                }
            }
        for(t in 2:T) {
            if(identical(dynamics, "autoreg"))
                N[i, t] <- N[i, t-1] * (omega[i, t-1] + gamma[i, t-1]) +
                    iota[i, t-1]
            else if(identical(dynamics, "trend"))
                N[i,t] <- N[i,t-1] * gamma[i,t-1] + iota[i, t-1]
            else if(identical(dynamics, "ricker"))
                N[i,t] <- N[i,t-1]*exp(gamma[i,t-1]*(1-N[i,t-1]/omega[i,t-1]))+
                    iota[i, t-1]
            else if(identical(dynamics, "gompertz"))
                N[i,1] <- N[i,t-1] * exp(gamma[i,t-1]*(1-log(N[i,t-1]+1)/
                  log(omega[i,t-1]+1))) + iota[i, t-1]
            else
                N[i,t] <- N[i,t-1] * omega[i,t-1] + gamma[i,t-1]
            if(delta[i, t] > 1) {
                for(d in 2:delta[i, t]) {
                    if(identical(dynamics, "autoreg"))
                        N[i, t] <- N[i, t] * (omega[i, t-1] + gamma[i, t-1]) +
                            iota[i, t-1]
                    else if(identical(dynamics, "trend"))
                        N[i, t] <- N[i, t] * gamma[i, t-1] + iota[i, t-1]
                    else if(identical(dynamics, "ricker"))
                        N[i, t] <- N[i, t] * exp(gamma[i, t-1] * (1 - N[i,t] /
                            omega[i,t-1]))+ iota[i, t-1]
                    else if(identical(dynamics, "gompertz"))
                        N[i, 1] <- N[i, t] * exp(gamma[i, t-1] * (1 -
                            log(N[i, t]+1) / log(omega[i, t-1] + 1))) +
                            iota[i, t-1]
                    else
                        N[i,t] <- N[i,t] * omega[i, t-1] + gamma[i, t-1]
                    }
                }
            }
        }
    N <- N[,rep(1:T, each=J)]


    p <- getP(object, na.rm)
    N * p
})


setMethod("fitted", "unmarkedFitOccuRN", function(object, K, na.rm = FALSE)
{
    data <- object@data
    des <- getDesign(data, object@formula, na.rm = na.rm)
    X <- des$X; V <- des$V
    X.offset <- des$X.offset; V.offset <- des$V.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    if (is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
        }
    y <- des$y	# getY(data) ... to be consistent w/NA handling?
    y <- truncateToBinary(y)
    M <- nrow(X)
    J <- ncol(y)
    lam <- exp(X %*% coef(object, 'state') + X.offset)
    r <- plogis(V %*% coef(object, 'det') + V.offset)
    if(missing(K))
        K <- object@K ##max(y, na.rm = TRUE) + 20
    lam <- rep(lam, each = J)
    fitted <- 1 - exp(-lam*r) ## analytical integration.
    return(matrix(fitted, M, J, byrow = TRUE))
})

setMethod("fitted", "unmarkedFitOccuMulti", function(object)
{

  S <- length(object@data@ylist)
  J <- ncol(object@data@ylist[[1]])
  pmat <- getP(object)

  fitted_list <- list()
  for (i in 1:S){
    marg_occ <- predict(object,'state',se.fit=F,species=i)$Predicted
    occmat <- t(tcrossprod(rep(1,J),marg_occ))
    fitted_list[[i]] <- pmat[[i]] * occmat
  }
  names(fitted_list) <- names(object@data@ylist)
  fitted_list

})

setMethod("fitted", "unmarkedFitOccuMS", function(object, na.rm = FALSE)
{
  data <- object@data
  T <- data@numPrimary
  J <- obsNum(data) / T
  N <- numSites(data)
  S <- data@numStates

  if(T>1){
    stop('Not implemented for dynamic models')
  }

  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=TRUE)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=TRUE)

  #Get predictions
  pr <- predict(object, 'psi', se.fit=F)
  pr <- sapply(pr,function(x) x$Predicted)
  pr <- pr[rep(1:nrow(pr),each=J),]

  pr_det <- predict(object, 'det', se.fit=F)
  pr_det <- sapply(pr_det,function(x) x$Predicted)

  fitvals <- rep(NA, nrow(pr_det))
  if(object@parameterization == 'multinomial'){
    pr <- cbind(1-rowSums(pr),pr)

    for (i in 1:nrow(pr_det)){
      occ <- pr[i,]
      sdp <- matrix(0,nrow=S,ncol=S)
      sdp[guide] <- pr_det[i,]
      sdp[,1] <- 1 - rowSums(sdp)
      fitvals[i] <- occ %*% sdp %*% 0:(S-1)
    }

  } else if(object@parameterization == 'condbinom'){
    stop('Conditional binomial parameterization not supported yet')
  }
  fit_out <- matrix(fitvals,N,J,byrow=T)

  fit_out
})

setMethod("fitted", "unmarkedFitColExt", function(object, na.rm = FALSE)
{
    data <- object@data
    M <- numSites(data)
    nY <- data@numPrimary
    J <- obsNum(data)/nY
    psiParms <- coef(object, 'psi')
    detParms <- coef(object, 'det')
    colParms <- coef(object, 'col')
    extParms <- coef(object, 'ext')
    formulaList <- list(psiformula=object@psiformula,
        gammaformula=object@gamformula,
        epsilonformula=object@epsformula,
        pformula=object@detformula)
    designMats <- getDesign(object@data, object@formula)
    V.itj <- designMats$V
    X.it.gam <- designMats$X.gam
    X.it.eps <- designMats$X.eps
    W.i <- designMats$W

    psiP <- plogis(W.i %*% psiParms)
    detP <- plogis(V.itj %*% detParms)
    colP <- plogis(X.it.gam  %*% colParms)
    extP <- plogis(X.it.eps %*% extParms)

    detP <- array(detP, c(J, nY, M))
    colP <- matrix(colP, M, nY, byrow = TRUE)
    extP <- matrix(extP, M, nY, byrow = TRUE)

    ## create transition matrices (phi^T)
    phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
    for(i in 1:M) {
        for(t in 1:(nY-1)) {
            phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t],
                1-extP[i,t]))
            }
        }

    ## first compute latent probs
    x <- array(NA, c(2, nY, M))
    x[1,1,] <- 1-psiP
    x[2,1,] <- psiP
    for(i in 1:M) {
        for(t in 2:nY) {
            x[,t,i] <- (phis[,,t-1,i] %*% x[,t-1,i])
            }
        }

    ## then compute obs probs
    fitted <- array(NA, c(J, nY, M))
    for(i in 1:M) {
        for(t in 1:nY) {
            for(j in 1:J) {
                fitted[j,t,i] <- (x[,t,i] %*%
                    matrix(c(1, 1 - detP[j,t,i], 0, detP[j,t,i]), 2, 2))[2]
                }
            }
        }

    return(matrix(fitted, M, J*nY, byrow = TRUE))
})



# This covers unmarkedFitGDS too
setMethod("fitted", "unmarkedFitGMM",
    function(object, na.rm = FALSE)
{

    # E[y_itj] = M_i * phi_it * cp_itj

    data <- object@data
    D <- getDesign(data, object@formula, na.rm = na.rm)
    Xlam <- D$Xlam
    Xphi <- D$Xphi
    Xdet <- D$Xdet

    Xlam.offset <- D$Xlam.offset
    Xphi.offset <- D$Xphi.offset
    Xdet.offset <- D$Xdet.offset
    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
    if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
    if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

    y <- D$y
    M <- nrow(y)
    T <- data@numPrimary
    J <- ncol(y) / T
    lambda <- drop(exp(Xlam %*% coef(object, 'lambda') + Xlam.offset))
    if(identical(object@mixture, "ZIP")) {
        psi <- plogis(coef(object, type="psi"))
        lambda <- (1-psi)*lambda
    }

    if(T==1)
        phi <- 1
    else
        phi <- plogis(Xphi %*% coef(object, 'phi') + Xphi.offset)
    phi.mat <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
    phi.ijt <- as.numeric(apply(phi.mat, 2, rep, times=J))
    cp <- getP(object, na.rm = na.rm)

    fitted <- lambda * phi.ijt * as.numeric(cp) # recycle
    fitted <- matrix(fitted, M, J*T)
    return(fitted)
})


setMethod("fitted", "unmarkedFitOccuTTD", function(object, na.rm = FALSE)
{

  N <- nrow(object@data@y)
  T <- object@data@numPrimary
  J <- ncol(object@data@y)/T

  #Get predicted values
  psi <- predict(object, 'psi', na.rm=FALSE)$Predicted
  psi <- cbind(1-psi, psi)
  est_p <- getP(object)
  est_p <- as.numeric(t(est_p))
  est_p <- cbind(1-est_p, est_p)

  if(T>1){
    p_col <- predict(object, 'col', na.rm=FALSE)$Predicted
    p_ext <- predict(object, 'ext', na.rm=FALSE)$Predicted
    rem_seq <- seq(T, length(p_col), T)
    p_col <- p_col[-rem_seq]
    p_ext <- p_ext[-rem_seq]
    phi <- cbind(1-p_col, p_col, p_ext, 1-p_ext)
  }

  ## first compute latent probs
  state <- array(NA, c(2, T, N))
  state[1:2,1,] <- t(psi)

  if(T>1){
    phi_ind <- 1
    for(n in 1:N) {
      for(t in 2:T) {
        phi_mat <- matrix(phi[phi_ind,], nrow=2, byrow=TRUE)
        state[,t,n] <- phi_mat %*% state[,t-1,n]
        phi_ind <- phi_ind + 1
      }
    }
  }

  ## then compute obs probs
  obs <- array(NA, c(J, T, N))
  p_ind <- 1
  for(n in 1:N) {
    for(t in 1:T) {
      for(j in 1:J) {
        pmat <- matrix(c(1,0, est_p[p_ind,]), nrow=2, byrow=TRUE)
        obs[j,t,n] <- (state[,t,n] %*% pmat)[2] #prob y=1
        p_ind <- p_ind + 1
      }
    }
  }

  matrix(obs, N, J*T, byrow=TRUE)
})


setMethod("fitted", "unmarkedFitNmixTTD", function(object, na.rm = FALSE)
{
  stop("This method is not implemented for nmixTTD at this time", call.=FALSE)
})

## # Identical to method for unmarkedFitGMM. Need to fix class structure
## setMethod("fitted", "unmarkedFitGPC",
##     function(object, na.rm = FALSE)
## {
##     data <- object@data
##     D <- getDesign(data, object@formula, na.rm = na.rm)
##     Xlam <- D$Xlam
##     Xphi <- D$Xphi
##     Xdet <- D$Xdet

##     Xlam.offset <- D$Xlam.offset
##     Xphi.offset <- D$Xphi.offset
##     Xdet.offset <- D$Xdet.offset
##     if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
##     if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
##     if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

##     y <- D$y
##     M <- nrow(y)
##     T <- data@numPrimary
##     J <- ncol(y) / T
##     lambda <- drop(exp(Xlam %*% coef(object, 'lambda') + Xlam.offset))
##     if(T==1)
##         phi <- 1
##     else
##         phi <- plogis(Xphi %*% coef(object, 'phi') + Xphi.offset)
##     phi.mat <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
##     phi.ijt <- as.numeric(apply(phi.mat, 2, rep, times=J))
##     cp <- getP(object, na.rm = na.rm)

##     fitted <- lambda * phi.ijt * as.numeric(cp) # recycle
##     fitted <- matrix(fitted, M, J*T)
##     return(fitted)
## })


setMethod("profile", "unmarkedFit",
    function(fitted, type, parm, seq)
{
    stopifnot(length(parm) == 1)
    MLE <- mle(fitted)
    SE(fitted[type])
    nPar <- length(mle(fitted))
    nll <- nllFun(fitted)

    ## create table to match parm numbers with est/parm numbers
    types <- names(fitted)
    numbertable <- data.frame(type = character(0), num = numeric(0))
    for(i in seq(length=length(types))) {
        length.est <- length(fitted[i]@estimates)
        numbertable <- rbind(numbertable,
        data.frame(type = rep(types[i], length.est),
            num = seq(length=length.est)))
        }
    parm.fullnums <- which(numbertable$type == type &
        numbertable$num == parm)

    f <- function(value) {
        fixedNLL <- genFixedNLL(nll, parm.fullnums, value)
        mleRestricted <- optim(rep(0,nPar), fixedNLL)$value
        mleRestricted
        }
        prof.out <- sapply(seq, f)
        prof.out <- cbind(seq, prof.out)
        new("profile", prof = prof.out)
})



setMethod("hessian", "unmarkedFit",
    function(object)
{
    object@opt$hessian
})


setMethod("update", "unmarkedFit",
    function(object, formula., ..., evaluate = TRUE)
{
    call <- object@call
    origFormula <- formula(call)
    if (is.null(call))
        stop("need an object with call slot")
    extras <- match.call(call=sys.call(-1),
                         expand.dots = FALSE)$...
    if (!missing(formula.)) {
        detformula <- as.formula(origFormula[[2]])
        stateformula <- as.formula(paste("~", origFormula[3], sep=""))
        newDetformula <- as.formula(formula.[[2]])
        upDetformula <- update.formula(detformula, newDetformula)
        newStateformula <- as.formula(paste("~", formula.[3], sep=""))
        upStateformula <- update.formula(stateformula, newStateformula)
        call$formula <- as.formula(paste(
			deparse(upDetformula, width.cutoff=500),
            deparse(upStateformula, width.cutoff=500)))
            }
    if (length(extras) > 0) {
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

setMethod("update", "unmarkedFitOccuMulti",
    function(object, detformulas, stateformulas, ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    if(!missing(detformulas)){
      call[["detformulas"]] <- detformulas
    } else {
      call[["detformulas"]] <- object@detformulas
    }
    if(!missing(stateformulas)){
      call[["stateformulas"]] <- stateformulas
    } else {
      call[["stateformulas"]] <- object@stateformulas
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

setMethod("update", "unmarkedFitOccuMS",
    function(object, detformulas, psiformulas, phiformulas, ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    if(!missing(detformulas)){
      call[["detformulas"]] <- detformulas
    } else {
      call[["detformulas"]] <- object@detformulas
    }
    if(!missing(psiformulas)){
      call[["psiformulas"]] <- psiformulas
    } else {
      call[["psiformulas"]] <- object@psiformulas
    }
    if(!missing(phiformulas)){
      call[["phiformulas"]] <- phiformulas
    } else {
      if(!is.null(call$phiformulas)){
        call[["phiformulas"]] <- object@phiformulas
      }
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

setMethod("update", "unmarkedFitOccuTTD",
    function(object, psiformula, gammaformula, epsilonformula,detformula,
             ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    if(!missing(psiformula)){
      call[["psiformula"]] <- psiformula
    }
    if(!missing(gammaformula)){
      call[["gammaformula"]] <- gammaformula
    }
    if(!missing(epsilonformula)){
      call[["epsilonformula"]] <- epsilonformula
    }
    if(!missing(detformula)){
      call[["detformula"]] <- detformula
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

setMethod("update", "unmarkedFitNmixTTD",
    function(object, stateformula ,detformula, ..., evaluate = TRUE)
{

    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    if(!missing(stateformula)){
      call[["stateformula"]] <- stateformula
    }
    if(!missing(detformula)){
      call[["detformula"]] <- detformula
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

setMethod("update", "unmarkedFitColExt",
    function(object, formula., ..., evaluate = TRUE)
{
    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
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



setMethod("update", "unmarkedFitGMM",
    function(object, lambdaformula, phiformula, pformula, ...,
        evaluate = TRUE)
{
    call <- object@call
    if (is.null(call))
        stop("need an object with call slot")
    formlist <- object@formlist
    if (!missing(lambdaformula))
        call$lambdaformula <- update.formula(formlist$lambdaformula,
            lambdaformula)
    if (!missing(phiformula))
        call$phiformula <- update.formula(formlist$phiformula,
            phiformula)
    if (!missing(pformula))
        call$pformula <- update.formula(formlist$pformula,
            pformula)
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


setMethod("update", "unmarkedFitDailMadsen",
    function(object, lambdaformula., gammaformula., omegaformula.,
        pformula., iotaformula., ..., evaluate = TRUE) {
    call <- object@call
    lambdaformula <- as.formula(call[['lambdaformula']])
    gammaformula <- as.formula(call[['gammaformula']])
    omegaformula <- as.formula(call[['omegaformula']])
    pformula <- as.formula(call[['pformula']])
    iotaformula <- as.formula(call[['iotaformula']])
    extras <- match.call(call=sys.call(-1),
                         expand.dots = FALSE)$...
    if (!missing(lambdaformula.)) {
        upLambdaformula <- update.formula(lambdaformula,
                                          lambdaformula.)
        call[['lambdaformula']] <- upLambdaformula
    }
    if (!missing(gammaformula.)) {
        upGammaformula <- update.formula(gammaformula, gammaformula.)
        call[['gammaformula']] <- upGammaformula
    }
    if (!missing(omegaformula.)) {
        upOmegaformula <- update.formula(omegaformula, omegaformula.)
        call[['omegaformula']] <- upOmegaformula
    }
    if (!missing(pformula.)) {
        upPformula <- update.formula(pformula, pformula.)
        call[['pformula']] <- upPformula
    }
    if (!missing(iotaformula.)) {
        upIotaformula <- update.formula(iotaformula, iotaformula.)
        call[['iotaformula']] <- upIotaformula
    }
    if (length(extras) > 0) {
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


setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))
setMethod("sampleSize", "unmarkedFit", function(object) {
    data <- getData(object)
    M <- numSites(data)
    M <- M - length(object@sitesRemoved)
    M
})


setGeneric("getData", function(object) standardGeneric("getData"))
setMethod("getData", "unmarkedFit",function(object) {
    object@data
})


setGeneric("nllFun", function(object) standardGeneric("nllFun"))
setMethod("nllFun", "unmarkedFit", function(object) object@nllFun)

setGeneric("mle", function(object) standardGeneric("mle"))
setMethod("mle", "unmarkedFit", function(object) object@opt$par)

setClass("profile", representation(prof = "matrix"))

setGeneric("smoothed",
    function(object, mean=TRUE) standardGeneric("smoothed"))
setMethod("smoothed","unmarkedFitColExt",
    function(object, mean) {
        if(mean) object@smoothed.mean
        else object@smoothed
    })

setGeneric("projected",
    function(object, mean=TRUE) standardGeneric("projected"))
setMethod("projected","unmarkedFitColExt", function(object, mean) {
    if(mean) object@projected.mean
    else object@projected
})

setMethod("plot", c("profile", "missing"), function(x) {
    plot(x@prof[,1], x@prof[,2], type = "l")
})


setMethod("residuals", "unmarkedFit", function(object, ...) {
    y <- getY(object@data)
    e <- fitted(object, na.rm = FALSE)
    r <- y - e
    return(r)
})

setMethod("residuals", "unmarkedFitOccu", function(object, ...) {
    y <- getY(object@data)
    y <- truncateToBinary(y)
    e <- fitted(object, na.rm = FALSE)
    r <- y - e
    return(r)
})

setMethod("residuals", "unmarkedFitOccuFP", function(object, ...) {
  cat("residuals is not implemented for occuFP at this time")
})


setMethod("residuals", "unmarkedFitOccuRN", function(object, ...) {
    y <- getY(object@data)
    y <- truncateToBinary(y)
    e <- fitted(object, na.rm = FALSE)
    r <- y - e
    return(r)
})

setMethod("residuals", "unmarkedFitOccuMulti", function(object, ...) {
  res_list <- list()
  ylist <- object@data@ylist
  fitlist <- fitted(object)

  for (i in seq_along(ylist)){
    res_list[[i]] <- ylist[[i]] - fitlist[[i]]
  }
  names(res_list) <- names(ylist)
  res_list
})

setMethod("residuals", "unmarkedFitOccuTTD", function(object, ...) {
  tmax <- object@data@surveyLength
  yraw <- object@data@y
  y <- ifelse(yraw<tmax,1,0)
  e <- fitted(object)
  y - e
})

setMethod("plot", c(x = "unmarkedFit", y = "missing"), function(x, y, ...)
{
    r <- residuals(x)
    e <- fitted(x, na.rm = FALSE)
    plot(e, r, ylab = "Residuals", xlab = "Predicted values")
    abline(h = 0, lty = 3, col = "gray")
})

setMethod("plot", c(x = "unmarkedFitOccuMulti", y = "missing"), function(x, y, ...)
{
  r <- do.call(rbind,residuals(x))
  e <- do.call(rbind,fitted(x))
  plot(e, r, ylab = "Residuals", xlab = "Predicted values")
  abline(h = 0, lty = 3, col = "gray")
})


setMethod("hist", "unmarkedFitDS", function(x, lwd=1, lty=1, ...) {
    ymat <- getY(getData(x))
    dbreaks <- getData(x)@dist.breaks
    nb <- length(dbreaks)
    mids <- (dbreaks[-1] - dbreaks[-nb]) / 2 + dbreaks[-nb]
    distances <- rep(mids, times=colSums(ymat))
    h <- hist(distances, plot=F, breaks=dbreaks)
    key <- x@keyfun
    survey <- x@data@survey
    switch(key,
           halfnorm = {
               sigma <- exp(coef(x, type="det"))
               if(length(sigma) > 1)
               stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- 2 * integrate(dnorm, dbreaks[1],
                                               dbreaks[nb], sd=sigma)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) 2 * dnorm(x, mean=0, sd=sigma),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drhn, dbreaks[1], dbreaks[nb],
                                           sigma=sigma)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drhn(r, sigma=sigma),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      })
           },
           exp = {		# This doesn't work on example fm4
               rate <- exp(coef(x, type="det"))
               if(length(rate) > 1)
                   stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- integrate(dxexp, dbreaks[1], dbreaks[nb],
                                           rate=rate)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) dxexp(x, rate=rate),
                               min(dbreaks),
                               max(dbreaks), add=T, lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drexp, dbreaks[1], dbreaks[nb],
                                           rate=rate)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drexp(r, rate=rate),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      })
           },
           hazard = {
               shape <- exp(coef(x, type="det"))
               scale <- exp(coef(x, type="scale"))
               if(length(shape) > 1)
                   stop("This method only works when there are no detection covars")
               switch(survey,
                      line = {
                          int <- integrate(dxhaz, dbreaks[1], dbreaks[nb],
                                           shape=shape, scale=scale)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(x) dxhaz(x, shape=shape,
                                                 scale=scale),
                               min(dbreaks), max(dbreaks), add=T,
                               lwd=lwd, lty=lty)
                      },
                      point = {
                          int <- integrate(drhaz, dbreaks[1], dbreaks[nb],
                                           shape=shape, scale=scale)$value
                          h$density <- h$density * int
                          plot(h, freq=F, ...)
                          plot(function(r) drhaz(r, shape=shape,
                                                 scale=scale),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      })
           },
           uniform = {
               switch(survey,
                      line = {
                          plot(h, freq=F, ...)
                          abline(h=1/max(dbreaks), lwd=lwd, lty=lty)
                      },
                      point = {
                          plot(h, freq=F, ...)
                          plot(function(r) (pi*r^2) / (pi*max(dbreaks)^2),
                               min(dbreaks), max(dbreaks), add=T, lwd=lwd,
                               lty=lty)
                      }
                      )}
           )
})




# ----------------------- CHILD CLASS METHODS ---------------------------

# Extract detection probs
setGeneric("getP", function(object, ...) standardGeneric("getP"))
setGeneric("getFP", function(object, ...) standardGeneric("getFP"))
setGeneric("getB", function(object, ...) standardGeneric("getB"))


setMethod("getP", "unmarkedFit", function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    V <- cbind(designMats$V, designMats$Z_det)
    V.offset <- designMats$V.offset
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))
    M <- nrow(y)
    J <- ncol(y)
    ppars <- coef(object, type = "det", fixedOnly=FALSE)
    p <- plogis(as.matrix(V %*% ppars + V.offset))
    p <- matrix(p, M, J, byrow = TRUE)
    return(p)
})


setMethod("getP", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  formula <- object@formula
  detformula <- object@detformula
  stateformula <- object@stateformula
  FPformula <- object@FPformula
  Bformula <- object@Bformula
  umf <- object@data
  designMats <- getDesign(umf, detformula,FPformula,Bformula,stateformula, na.rm = na.rm)
  y <- designMats$y
  V <- designMats$V
  V.offset <- designMats$V.offset
  if (is.null(V.offset))
    V.offset <- rep(0, nrow(V))
  M <- nrow(y)
  J <- ncol(y)
  ppars <- coef(object, type = "det")
  p <- plogis(V %*% ppars + V.offset)
  p <- matrix(p, M, J, byrow = TRUE)
  return(p)
})

setMethod("getP", "unmarkedFitOccuMulti", function(object)
{

  ylist <- object@data@ylist
  S <- length(ylist)
  N <- nrow(ylist[[1]])
  maxOrder <- object@call$maxOrder
  if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)
  dm <- getDesign(object@data,object@detformulas,object@stateformulas, maxOrder=maxOrder)
  pred <- predict(object,'det',se.fit=F)
  dets <- do.call(cbind,lapply(pred,`[`,,1))
  #ugly mess
  out <- list()
  for (i in 1:S){
    pmat <- array(NA,dim(ylist[[1]]))
    for (j in 1:N){
      ps <- dets[dm$yStart[j]:dm$yStop[j],i]
      not_na <- !is.na(ylist[[i]][j,])
      pmat[j,not_na] <- ps
    }
    out[[i]] <- pmat
  }
  names(out) <- names(ylist)
  out
})

setMethod("getP", "unmarkedFitOccuMS", function(object)
{
  J <- ncol(object@data@y)
  N <- nrow(object@data@y)
  pred <- predict(object, 'det', se.fit=F)
  lapply(pred, function(x) matrix(x$Predicted, nrow=N, ncol=J, byrow=T))
})

setMethod("getP", "unmarkedFitOccuTTD", function(object)
{

  N <- nrow(object@data@y)
  lam <- predict(object, 'det', na.rm=FALSE)$Predicted
  tmax <- as.numeric(t(object@data@surveyLength))
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")

  not_na <- !is.na(lam)
  est_p <- rep(NA, length(lam))
  if(tdist == "weibull"){
    k <- exp(coef(object)['k(k)'])
    est_p[not_na] <- stats::pweibull(tmax[not_na], k, 1/lam[not_na])
  } else {
    est_p[not_na] <- stats::pexp(tmax[not_na], lam[not_na])
  }

  matrix(est_p, nrow=N, byrow=TRUE)
})

setMethod("getFP", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  formula <- object@formula
  detformula <- object@detformula
  stateformula <- object@stateformula
  FPformula <- object@FPformula
  Bformula <- object@Bformula
  umf <- object@data
  designMats <- getDesign(umf, detformula,FPformula,Bformula,stateformula, na.rm = na.rm)
  type = object@type
  y <- designMats$y
  U <- designMats$U
  U.offset <- designMats$U.offset
  if (is.null(U.offset))
    U.offset <- rep(0, nrow(U))
  M <- nrow(y)
  J <- ncol(y)
  fpars <- coef(object, type = "fp")
  f <- plogis(U %*% fpars + U.offset)
  f <- matrix(f, M, J, byrow = TRUE)
  if (type[1]!=0){
    f[,1:type[1]] = 0
  }
  return(f)
})

setMethod("getB", "unmarkedFitOccuFP", function(object, na.rm = TRUE)
{
  formula <- object@formula
  detformula <- object@detformula
  stateformula <- object@stateformula
  FPformula <- object@FPformula
  Bformula <- object@Bformula
  umf <- object@data
  designMats <- getDesign(umf, detformula,FPformula,Bformula,stateformula, na.rm = na.rm)
  y <- designMats$y
  W <- designMats$W
  W.offset <- designMats$W.offset
  if (is.null(W.offset))
    W.offset <- rep(0, nrow(W))
  M <- nrow(y)
  J <- ncol(y)
  type = object@type
  if (type[3]!=0){
    bpars <- coef(object, type = "b")
  b <- plogis(W %*% bpars + W.offset)
  b <- matrix(b, M, J, byrow = TRUE)
  }
  if (type[3]==0){
    b <- matrix(0, M, J)
  }
  return(b)
})


setMethod("getP", "unmarkedFitDS",
    function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    V <- cbind(designMats$V, designMats$Z_det)
    V.offset <- designMats$V.offset
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))
    M <- nrow(y)
    J <- ncol(y)
    ppars <- coef(object, type = "det", fixedOnly=FALSE)
    db <- umf@dist.breaks
    w <- diff(db)
    survey <- umf@survey
    key <- object@keyfun
    tlength <- umf@tlength

    cp <- u <- a <- matrix(NA, M, J)
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
            if(J > 1){
              for(j in 2:J)
                a[i, j] <- pi*db[j+1]^2 - sum(a[i, 1:(j-1)])
            }
            u[i,] <- a[i,] / sum(a[i,])
        }
            })


    switch(key,
        halfnorm = {
            sigma <- exp(V %*% ppars + V.offset)
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
                        cp[i, j] <- integrate(grhn, db[j], db[j+1],
                            sigma=sigma[i], rel.tol=1e-4)$value *
                            2 * pi / a[i, j]
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            },
        exp = {
            rate <- exp(V %*% ppars + V.offset)
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        cp[i, j] <- integrate(gxexp, db[j], db[j+1],
                            rate=rate[i], rel.tol=1e-4)$value / w[j]
                        }},
                point = {
                    for(j in 1:J) {
                        cp[i, j] <- integrate(grexp, db[j], db[j+1],
                            rate=rate[i], rel.tol=1e-4)$value *
                            2 * pi / a[i, j]
                        }
                    })
                cp[i,] <- cp[i,] * u[i,]
                }
            },
        hazard = {
            shape <- exp(V %*% ppars + V.offset)
            scale <- exp(coef(object, type="scale"))
            for(i in 1:M) {
                switch(survey,
                line = {
                    for(j in 1:J) {
                        cp[i, j] <- integrate(gxhaz, db[j], db[j+1],
                            shape=shape[i], scale=scale,
                            rel.tol=1e-4)$value / w[j]
                        }},
                point = {
                    for(j in 1:J) {
                        cp[i, j] <- integrate(grhaz, db[j], db[j+1],
                            shape = shape[i], scale=scale,
                            rel.tol=1e-4)$value * 2 * pi / a[i, j]
                    }})
                cp[i,] <- cp[i,] * u[i,]
                }
            },
		uniform = cp <- u)
    return(cp)
})




# Should this return p or pi. Right now it's pi without phi.
setMethod("getP", "unmarkedFitGDS",
    function(object, na.rm = TRUE)
{
#    browser()
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    Xdet <- designMats$Xdet
    Xdet.offset <- designMats$Xdet.offset
    if (is.null(Xdet.offset))
        Xdet.offset <- rep(0, nrow(Xdet))
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    ppars <- coef(object, type = "det")
    db <- umf@dist.breaks
    w <- diff(db)
    survey <- umf@survey
    key <- object@keyfun
    tlength <- umf@tlength

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


    cp <- array(NA, c(M, J, T))
    switch(key,
        halfnorm = {
            sigma <- exp(Xdet %*% ppars + Xdet.offset)
            sigma <- matrix(sigma, M, T, byrow=TRUE)
            for(i in 1:M) {
                for(t in 1:T) {
                    if(is.na(sigma[i,t])) next
                    switch(survey,
                    line = {
                        f.0 <- 2 * dnorm(0, 0, sd=sigma[i, t])
                        int <- 2 * (pnorm(db[-1], 0, sd=sigma[i, t]) -
                            pnorm(db[-(J+1)], 0, sd=sigma[i, t]))
                        cp[i,,t] <- int / f.0 / w
                        },
                    point = {
                        for(j in 1:J) {
                            cp[i, j, t] <- integrate(grhn, db[j], db[j+1],
                                sigma=sigma[i, t], rel.tol=1e-4)$value *
                                2 * pi / a[i, j]
                            }
                        })
                    cp[i,,t] <- cp[i,,t] * u[i,]
                    }
                }
            },
        exp = {
            rate <- exp(Xdet %*% ppars + Xdet.offset)
            rate <- matrix(rate, M, T, byrow=TRUE)
            for(i in 1:M) {
                for(t in 1:T) {
                if(is.na(rate[i,t])) next
                switch(survey,
                line = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(gxexp, db[j], db[j+1],
                            rate=rate[i,t], rel.tol=1e-4)$value / w[j]
                        }},
                point = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(grexp, db[j], db[j+1],
                            rate=rate[i,t], rel.tol=1e-4)$value *
                            2 * pi / a[i, j]
                        }
                    })
                cp[i,,t] <- cp[i,,t] * u[i,]
                }}
            },
        hazard = {
            shape <- exp(Xdet %*% ppars + Xdet.offset)
            shape <- matrix(shape, M, T, byrow=TRUE)
            scale <- exp(coef(object, type="scale"))
            for(i in 1:M) {
                for(t in 1:T) {
                if(is.na(shape[i,t])|is.na(scale)) next
                switch(survey,
                line = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(gxhaz, db[j], db[j+1],
                            shape=shape[i,t], scale=scale,
                            rel.tol=1e-4)$value / w[j]
                        }},
                point = {
                    for(j in 1:J) {
                        cp[i, j, t] <- integrate(grhaz, db[j], db[j+1],
                            shape = shape[i,t], scale=scale,
                            rel.tol=1e-4)$value * 2 * pi / a[i, j]
                    }})
                cp[i,,t] <- cp[i,,t] * u[i,]
                }}
            },
	uniform = {
#            browser()
            cp[] <- u
        })
    cp <- matrix(cp, nrow=M)
    return(cp)
})


setMethod("getP", "unmarkedFitDSO",
    function(object, na.rm = TRUE)
{
    umf <- getData(object)
    y <- getY(umf)
    M <- numSites(umf)

    keep <- 1:M
    if(length(object@sitesRemoved)>0){
      keep <- keep[-object@sitesRemoved]
    }

    M <- length(keep)
    T <- umf@numPrimary
    J <- ncol(y) / T

    sig <- matrix(NA, M, T)
    if(object@keyfun != "uniform"){
      sig <- predict(object, type="det")$Predicted
      sig <- matrix(sig, M, T, byrow=TRUE)
    }

    scale <- 0.0
    if(object@keyfun == "hazard"){
      scale <- backTransform(object, type="scale")@estimate
    }

    db <- umf@dist.breaks
    w <- diff(db)
    ua <- getUA(umf)
    u <- ua$u[keep,]; a <- ua$a[keep,]

    cp <- array(NA, c(M, J, T))
    for (i in 1:M){
      for (t in 1:T){
        cp[i,,t] <- getDistCP(object@keyfun, sig[i,t], scale, umf@survey,
                              db, w, a[i,], u[i,])
      }
    }
    matrix(cp, nrow=M)
})


setMethod("getP", "unmarkedFitMPois", function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    piFun <- object@data@piFun
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    V <- as.matrix(cbind(designMats$V, designMats$Z_det))
    V.offset <- designMats$V.offset
    if (is.null(V.offset))
        V.offset <- rep(0, nrow(V))
    M <- nrow(y)
    J <- obsNum(umf) #ncol(y)
    ppars <- coef(object, type = "det", fixedOnly=FALSE)
    p <- plogis(V %*% ppars + V.offset)
    p <- matrix(p, M, J, byrow = TRUE)
    pi <- do.call(piFun, list(p = p))
    return(pi)
})


setMethod("getP", "unmarkedFitMMO", function(object, na.rm = TRUE)
{

  umf <- object@data

  D <- getDesign(umf, object@formula, na.rm=na.rm)
  beta <- coef(object, type='det')
  off <- D$Xp.offset
  if(is.null(off)) off <- rep(0, nrow(D$Xp))
  plong <- plogis(D$Xp %*% beta + off)

  M <- nrow(D$y)
  T <- umf@numPrimary
  J <- ncol(getY(umf)) / T

  pmat <- aperm(array(plong, c(J,T,M)), c(3,1,2))

  pout <- array(NA, c(M,J,T))
  for (t in 1:T){
    pout[,,t] <- do.call(umf@piFun, list(p=pmat[,,t]))
  }
  matrix(aperm(pout,c(2,3,1)), M, J*T, byrow=TRUE)

})


setMethod("getP", "unmarkedFitPCO", function(object, na.rm = TRUE)
{
    umf <- object@data
    D <- getDesign(umf, object@formula, na.rm = na.rm)
    y <- D$y
    Xp <- D$Xp
    Xp.offset <- D$Xp.offset
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T*J)
    ppars <- coef(object, type = "det")
    p <- plogis(Xp %*% ppars + Xp.offset)
    p <- matrix(p, M, J*T, byrow = TRUE)
    return(p)
})



setMethod("getP", "unmarkedFitColExt", function(object, na.rm = TRUE)
{
    data <- object@data
    detParms <- coef(object, 'det')
    D <- getDesign(object@data, object@formula, na.rm=na.rm)
    y <- D$y
    V <- D$V

    M <- nrow(y)	# M <- nrow(X.it)
    nY <- data@numPrimary
    J <- obsNum(data)/nY

    p <- plogis(V %*% detParms)
    p <- array(p, c(J, nY, M))
    p <- aperm(p, c(3, 1, 2))
    p <- matrix(p, nrow=M)
    return(p)
})



setMethod("getP", "unmarkedFitGMM",
    function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- object@formlist$pformula
    piFun <- object@data@piFun
    umf <- object@data
    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xdet <- D$Xdet
    Xdet.offset <- D$Xdet.offset
    if (is.null(Xdet.offset))
        Xdet.offset <- rep(0, nrow(Xdet))

    M <- nrow(y)
    T <- object@data@numPrimary
    R <- numY(object@data) / T
    J <- obsNum(object@data) / T

    ppars <- coef(object, type = "det")
    p <- plogis(Xdet %*% ppars + Xdet.offset)
    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, J, T))
    p <- aperm(p, c(1,3,2))

    cp <- array(as.numeric(NA), c(M, T, R))
    for(t in 1:T) cp[,t,] <- do.call(piFun, list(matrix(p[,t,], M, J)))
    cp <- aperm(cp, c(1,3,2))
    cp <- matrix(cp, nrow=M, ncol=numY(object@data))

    return(cp)
})


setMethod("getP", "unmarkedFitGPC",
    function(object, na.rm = TRUE)
{
    formula <- object@formula
    detformula <- object@formlist$pformula
    umf <- object@data
    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xdet <- D$Xdet
    Xdet.offset <- D$Xdet.offset
    if (is.null(Xdet.offset))
        Xdet.offset <- rep(0, nrow(Xdet))
    R <- nrow(y)
    T <- object@data@numPrimary
    J <- ncol(y) / T
    ppars <- coef(object, type = "det")
    p <- plogis(Xdet %*% ppars + Xdet.offset)
    p <- matrix(p, nrow=R, byrow=TRUE)
    return(p)
})

#Y extractors for unmarkedFit objects
setMethod("getY", "unmarkedFit", function(object) object@data@y)
setMethod("getY", "unmarkedFitOccu", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY", "unmarkedFitOccuRN", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY", "unmarkedFitColExt", function(object) {
            truncateToBinary(object@data@y)
})
setMethod("getY", "unmarkedFitOccuMulti", function(object) {
            object@data@ylist
})



setMethod("simulate", "unmarkedFitDS",
    function(object, nsim = 1, seed = NULL, na.rm=TRUE)
{
    formula <- object@formula
    umf <- object@data
    db <- umf@dist.breaks
    w <- diff(db)
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    X <- as.matrix(cbind(designMats$X, designMats$Z_state))
    X.offset <- designMats$X.offset
    if (is.null(X.offset))
        X.offset <- rep(0, nrow(X))
    M <- nrow(y)
    J <- ncol(y)
    lamParms <- coef(object, type = "state", fixedOnly=FALSE)
    lambda <- drop(exp(X %*% lamParms + X.offset))
    if(identical(object@output, "density")) {
        switch(umf@survey,
            line = {
                tlength <- umf@tlength
                A <- tlength * max(db) * 2
                },
            point = {
                A <- pi * max(db)^2
                })
        switch(umf@unitsIn,
            m = A <- A / 1e6,
            km = A <- A)
        switch(object@unitsOut,
            m = A <- A * 1e6,
            ha = A <- A * 100,
            kmsq = A <- A)
        lambda <- lambda * A
        }
    cp <- getP(object, na.rm = na.rm)
    simList <- list()
    for(i in 1:nsim) {
        yvec <- rpois(M * J, lambda * cp)
        simList[[i]] <- matrix(yvec, M, J)
        }
    return(simList)
})



setMethod("simulate", "unmarkedFitPCount",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    formula <- object@formula
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    X <- cbind(designMats$X, designMats$Z_state)
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    lamParms <- coef(object, type = "state", fixedOnly=FALSE)
    lam <- as.numeric(exp(X %*% lamParms + X.offset))
    lamvec <- rep(lam, each = J)
    pvec <- c(t(getP(object, na.rm = na.rm)))
    mix <- object@mixture
    simList <- list()
    for(i in 1:nsim) {
        switch(mix,
               P = N <- rpois(M, lam),
               NB = N <- rnbinom(M, size = exp(coef(object["alpha"])),
                                 mu = lam),
               ZIP = {
                   psi <- plogis(coef(object["psi"]))
                   N <- rzip(M, lam, psi)
               }
            )
        yvec <- rbinom(M * J, size = rep(N, each = J), prob = pvec)
        simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
        }
    return(simList)
})


#Simulate open-population abundance
simOpenN <- function(object, na.rm)
{
    mix <- object@mixture
    dynamics <- object@dynamics
    umf <- object@data
    #To partially handle old saved model objects
    fix <- tryCatch(object@fix, error=function(e) "none")
    immigration <- tryCatch(object@immigration, error=function(e) FALSE)
    D <- getDesign(umf, object@formula, na.rm = na.rm)
    Xlam <- D$Xlam; Xgam <- D$Xgam; Xom <- D$Xom; Xp <- D$Xp; Xiota <- D$Xiota
    Xlam.offset <- D$Xlam.offset; Xgam.offset <- D$Xgam.offset
    Xom.offset <- D$Xom.offset; Xp.offset <- D$Xp.offset
    Xiota.offset <- D$Xiota.offset
    delta <- D$delta

    y <- D$y
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T

    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
    if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
    if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
    if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T*J)
    if(is.null(Xiota.offset)) Xiota.offset <- rep(0, M*(T-1))

    lambda <- drop(exp(Xlam %*% coef(object, 'lambda') + Xlam.offset))
    if(fix == "gamma"){
        gamma <- matrix(0, M, T-1)
    } else if(dynamics != "notrend"){
        gamma <- matrix(exp(Xgam %*% coef(object, 'gamma') + Xgam.offset),
                        M, T-1, byrow=TRUE)
    } else {
        gamma <- matrix(NA, M, T-1)
    }
    if (fix == "omega")
        omega <- matrix(1, M, T-1)
    else if(dynamics == "trend")
        omega <- matrix(-999, M, T-1) # placeholder
    else if(identical(dynamics, "ricker"))
        omega <- matrix(exp(Xom %*% coef(object, 'omega') + Xom.offset),
                        M, T-1, byrow=TRUE)
    else if(identical(dynamics, "gompertz"))
        omega <- matrix(exp(Xom %*% coef(object, 'omega') + Xom.offset),
                        M, T-1, byrow=TRUE)
    else
        omega <- matrix(plogis(Xom %*% coef(object, 'omega') + Xom.offset),
                        M, T-1, byrow=TRUE)
    if(immigration)
        iota <- matrix(exp(Xiota %*% coef(object, 'iota') + Xiota.offset),
                        M, T-1, byrow=TRUE)
    else
        iota <- matrix(0, M, T-1)
    if(identical(mix, "ZIP"))
        psi <- plogis(coef(object, type="psi"))

    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    for(i in 1:M) {
            switch(mix,
                   P = N[i, 1] <- rpois(1, lambda),
                   NB = N[i, 1] <- rnbinom(1, size =
                       exp(coef(object["alpha"])), mu = lambda),
                   ZIP = N[i,1] <- rzip(1, lambda[i], psi))
            if(delta[i, 1] > 1) {
                for(d in 2:delta[i, 1]) {
                    if(dynamics == "trend")
                        N[i,1] <- rpois(1, N[i,1]*gamma[i,1]+iota[i,1])
                    else if(dynamics == "ricker")
                        N[i,1] <- rpois(1, N[i,1]*exp(gamma[i, 1] * (1 -
                          N[i, 1] / omega[i, 1])) + iota[i, 1])
                    else if(dynamics == "gompertz")
                        N[i,1] <- rpois(1, N[i, 1] * exp(gamma[i, 1] * (1 -
                          log(N[i, 1] + 1)/log(omega[i, 1] + 1))) +
                          iota[i, 1])
                    else if(dynamics == "autoreg")
                        N[i,1] <- rbinom(1, N[i,1], omega[i,1]) +
                            rpois(1, gamma[i,1]*N[i,1] + iota[i, 1])
                    else
                        N[i,1] <- rbinom(1, N[i,1], omega[i,1]) +
                            rpois(1, gamma[i, 1])
                }
            }
            # Might need more NA handling here...ignore gaps?
            for(t in 2:T) {
                if(is.na(omega[i, t-1]) | is.na(gamma[i,t-1]))
                    N[i, t] <- N[i, t-1] # just a place holder
                else{
                    if(!identical(dynamics, "trend"))
                        S[i, t-1] <- rbinom(1, N[i, t-1], omega[i, t-1])
                    if(identical(dynamics, "autoreg"))
                        gamma[i, t-1] <- gamma[i, t-1] * N[i,t-1] + iota[i,t-1]
                    else if(identical(dynamics, "notrend"))
                        gamma[i, t-1] <- (1-omega[i, t-1]) * lambda[i]
                    G[i, t-1] <- rpois(1, gamma[i, t-1])
                    if(identical(dynamics, "trend"))
                        N[i,t] <- rpois(1, N[i,t-1]*gamma[i,t-1]+iota[i,t-1])
                    else if(identical(dynamics, "ricker"))
                        N[i,t] <- rpois(1, N[i,t-1]*exp(gamma[i,t-1]*(1-N[i,
                            t-1]/omega[i,t-1]))+iota[i,t-1])
                    else if(identical(dynamics, "gompertz"))
                        N[i,t] <- rpois(1, N[i,t-1]*exp(gamma[i,t-1]*(1-
                            log(N[i,t-1] + 1) / log(omega[i,t-1] + 1))) +
                            iota[i,t-1])
                    else
                        N[i, t] <- S[i, t-1] + G[i, t-1]
                    if(delta[i, t] > 1) {
                        for(d in 2:delta[i, 1]) {
                            if(dynamics == "trend")
                                N[i,t] <- rpois(1, N[i,t]*gamma[i,t-1]+
                                    iota[i,t-1])
                            else if(identical(dynamics, "ricker"))
                                N[i,t] <- rpois(1, N[i,t]*exp(gamma[i,t-1]*(1-
                                    N[i,t]/omega[i,t-1]))+iota[i,t-1])
                            else if(identical(dynamics, "gompertz"))
                                N[i,t] <- rpois(1, N[i,t]*exp(gamma[i,t-1]*(1-
                                    log(N[i,t] + 1)/
                                    log(omega[i,t-1] + 1))) + iota[i,t-1])
                            else {
                                S[i,t-1] <- rbinom(1, N[i,t], omega[i,t-1])
                                G[i,t-1] <- rpois(1, gamma[i, t-1])
                                N[i, t] <- S[i, t-1] + G[i, t-1]
                            }
                        }
                    }
                }
            }
        }
    N
}

setMethod("simulate", "unmarkedFitPCO",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{

    umf <- object@data
    M <- numSites(umf)
    T <- umf@numPrimary
    J <- ncol(getY(umf)) / T
    D <- getDesign(umf, object@formula, na.rm = na.rm)
    y <- D$y
    y.na <- is.na(y)
    p <- getP(object, na.rm = na.rm)
    simList <- list()
    for(s in 1:nsim) {
        y.sim <- matrix(NA, M, J*T)
        N <- simOpenN(object, na.rm)
        N <- N[,rep(1:T, each=J)]
        y.sim[!y.na] <- rbinom(sum(!y.na), N[!y.na], p[!y.na])
        simList[[s]] <- y.sim
    }
    return(simList)
})


setMethod("simulate", "unmarkedFitDailMadsen",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
  umf <- object@data
  D <- getDesign(umf, object@formula, na.rm = na.rm)
  y <- D$y
  y.na <- is.na(y)
  M <- numSites(umf)
  T <- umf@numPrimary
  J <- ncol(getY(umf)) / T
  simList <- list()
  p <- getP(object, na.rm = na.rm)
  p <- array(p, c(M,J,T))
  cp <- array(NA, c(M,J+1,T))
  for (i in 1:M){
    for (t in 1:T){
      cp[i, 1:J, t] <- p[i,,t]
      cp[i, J+1, t] <- 1 - sum(p[i,,t], na.rm=TRUE)
    }
  }

  for(s in 1:nsim) {
    y.sim <- matrix(NA, M, J*T)
    N <- simOpenN(object, na.rm)

    for(i in 1:M) {
      yst <- 1
        for(t in 1:T) {
          yend <- yst + J - 1
          #rmultinom2 in utils.R
          y.it <- as.integer(rmultinom2(1, N[i,t], prob=cp[i,,t]))
          y.sim[i,yst:yend] <- y.it[1:J]
          yst <- yst + J
        }
    }
    y.sim[y.na] <- NA
    simList[[s]] <- y.sim
  }
  return(simList)
})


setMethod("simulate", "unmarkedFitMPois",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    formula <- object@formula
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    X <- as.matrix(cbind(designMats$X, designMats$Z_state))
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    lamParms <- coef(object, type = "state", fixedOnly=FALSE)
    lam <- as.numeric(exp(X %*% lamParms + X.offset))
    lamvec <- rep(lam, each = J)
    pivec <- as.vector(t(getP(object, na.rm = na.rm)))
    simList <- list()
    for(i in 1:nsim) {
        yvec <- rpois(M * J, lamvec * pivec)
        simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
        }
    return(simList)
})




setMethod("simulate", "unmarkedFitOccu",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    formula <- object@formula
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    X <- cbind(designMats$X, designMats$Z_state)
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    psiParms <- coef(object, type = "state", fixedOnly=FALSE)
    psi <- as.numeric(plogis(as.matrix(X %*% psiParms + X.offset)))
    p <- c(t(getP(object,na.rm = na.rm)))
    simList <- list()
    for(i in 1:nsim) {
        Z <- rbinom(M, 1, psi)
        Z[object@knownOcc] <- 1
        yvec <- rep(Z, each = J)*rbinom(M * J, 1, prob = p)
            simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
        }
    return(simList)
})


setMethod("simulate", "unmarkedFitOccuFP",
          function(object, nsim = 1, seed = NULL, na.rm = TRUE)
          {
            detformula <- object@detformula
            stateformula <- object@stateformula
            FPformula <- object@FPformula
            Bformula <- object@Bformula
            umf <- object@data
            designMats <- getDesign(umf, detformula,FPformula,Bformula,stateformula, na.rm = na.rm)
            X <- designMats$X; V <- designMats$V; U <- designMats$U; W <- designMats$W;
            y <- designMats$y
            X.offset <- designMats$X.offset; V.offset <- designMats$V.offset; U.offset <- designMats$U.offset; W.offset <- designMats$W.offset
            if(is.null(X.offset)) {
              X.offset <- rep(0, nrow(X))
            }
            if(is.null(V.offset)) {
              V.offset <- rep(0, nrow(V))
            }
            if(is.null(U.offset)) {
              U.offset <- rep(0, nrow(U))
            }
            if(is.null(W.offset)) {
              W.offset <- rep(0, nrow(W))
            }
            M <- nrow(y)
            J <- ncol(y)
            allParms <- coef(object, altNames = FALSE)
            psiParms <- coef(object, type = "state")
            psi <- as.numeric(plogis(X %*% psiParms + X.offset))
            p <- c(t(getP(object)))
            fp <- c(t(getFP(object)))
            b <- c(t(getB(object)))
            simList <- list()
            for(i in 1:nsim) {
              Z <- rbinom(M, 1, psi)
              Z[object@knownOcc] <- 1
              Z <- rep(Z, each = J)
              P <- matrix(0,M*J,3)
              P[,1] <- Z*rbinom(M * J, 1, prob = (1-p)) + (1-Z)*rbinom(M * J, 1, prob = (1-fp))
              P[,2] <- (1-P[,1])*(1-Z) + (1-P[,1])*rbinom(M * J, 1, prob = (1-b))*Z
              P[,3] <- 1 - P[,1]-P[,2]
              yvec <- sapply(1:(M*J),function(x) which(as.logical(rmultinom2(1,1,P[x,])))-1)
              simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
            }
            return(simList)
          })


setMethod("simulate", "unmarkedFitOccuMulti",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    data <- object@data
    ynames <- names(object@data@ylist)
    maxOrder <- object@call$maxOrder
    if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)
    dm <- getDesign(object@data, object@detformulas, object@stateformulas, maxOrder)
    psi <- predict(object, "state",se.fit=F)$Predicted
    p <- getP(object)

    simList <- list()
    for (s in 1:nsim){
      #True state
      ztruth <- matrix(NA,nrow=dm$N,ncol=dm$S)
      nZ <- nrow(dm$z)
      for (i in 1:dm$N){
        ztruth[i,] <- as.matrix(dm$z[sample(nZ,1,prob=psi[i,]),])
      }

      y <- list()
      for (i in 1:dm$S){
        y[[i]] <- matrix(NA,dm$N,dm$J)
        for (j in 1:dm$N){
          for (k in 1:dm$J){
            if(!is.na(p[[i]][j,k])){
              y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p[[i]][j,k])
            }
          }
        }
      }
      names(y) <- ynames
      simList[[s]] <- y
    }
    simList
})


setMethod("simulate", "unmarkedFitOccuMS",
    function(object, nsim = 1, seed = NULL, na.rm=TRUE)
{

  S <- object@data@numStates
  N <- numSites(object@data)
  T <- object@data@numPrimary
  J <- obsNum(object@data) / T

  prm <- object@parameterization
  psi_raw <- predict(object, "psi", se.fit=F)
  psi_raw <- sapply(psi_raw, function(x) x$Predicted)
  p <- getP(object)

  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=TRUE)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=TRUE)

  out <- vector("list",nsim)

  for (i in 1:nsim){

  #State process
  if(prm == "multinomial"){
    psi <- cbind(1-apply(psi_raw,1,sum),psi_raw)
  } else if (prm == "condbinom"){
    psi <- matrix(NA, nrow=N, ncol=S)
    psi[,1] <- 1-psi_raw[,1]
    psi[,2] <- (1-psi_raw[,2])*psi_raw[,1]
    psi[,3] <- psi_raw[,1]*psi_raw[,2]
  }

  z <- matrix(NA, nrow=N, ncol=T)

  #initial occupancy
  for (n in 1:N){
    z[n,1] <- sample(0:(S-1), 1, prob=psi[n,])
  }

  #transitions if T>1----------------------------------------------------------
  get_phimat <- function(prob_vec){
    if(prm=="multinomial"){
      out <- matrix(NA, nrow=S, ncol=S)
      out[outer(1:S, 1:S, function(i,j) i!=j)] <- prob_vec
      out <- t(out)
      diag(out) <- 1 - rowSums(out,na.rm=T)
      return(out)
    } else if(prm == "condbinom"){
      out <- matrix(prob_vec, nrow=S)
      return(cbind(1-out[,1], out[,1]*(1-out[,2]), out[,1]*out[,2]))
    }
  }

  if(T>1){
    phi_raw <- predict(object, "phi", se.fit=F)
    phi_raw <- sapply(phi_raw, function(x) x$Predicted)
    phi_index <- 1
    for (n in 1:N){
      for (t in 2:T){
        phimat <- get_phimat(phi_raw[phi_index,])
        phi_t <- phimat[(z[n,(t-1)]+1),]
        z[n,t] <- sample(0:(S-1), 1, prob=phi_t)
        phi_index <- phi_index+1
      }
    }
  }
  #----------------------------------------------------------------------------

  #Detection process
  y <- matrix(0, nrow=N, ncol=J*T)
  for (n in 1:N){
    yindex <- 1
    for (t in 1:T){
      for (j in 1:J){
        if(prm == "multinomial"){
          probs_raw <- sapply(p, function(x) x[n,yindex])
          # Make sure output is NA if probs have NA
          if(any(is.na(probs_raw))){
            y[n,yindex] <- NA
            yindex <- yindex + 1
            next
          }

          sdp <- matrix(0, nrow=S, ncol=S)
          sdp[guide] <- probs_raw
          sdp[,1] <- 1 - rowSums(sdp)

          probs <- sdp[z[n,t]+1,]

        } else if (prm == "condbinom"){
          p11 <- p[[1]][n,yindex]
          p12 <- p[[2]][n,yindex]
          p22 <- p[[3]][n,yindex]
          # Trap NAs in probability of detection
          if(any(is.na(c(p11, p12, p22)))){
            y[n,yindex] <- NA
            next
          }
          probs <- switch(z[n,t]+1,
                          c(1,0,0),
                          c(1-p11,p11,0),
                          c(1-p12,p12*(1-p22),p12*p22))
        }
        # this NA trap probably isn't necessary but leaving it in just in case
        if(all(!is.na(probs))){
          y[n,yindex] <- sample(0:(S-1), 1, prob=probs)
        } else {
          y[n,yindex] <- NA
        }
        yindex <- yindex + 1
      }
    }
  }

  out[[i]] <- y
  }

  out
})


setMethod("simulate", "unmarkedFitOccuTTD",
          function(object,  nsim = 1, seed = NULL, na.rm = FALSE)
{

  N <- nrow(object@data@y)
  T <- object@data@numPrimary
  J <- ncol(object@data@y)/T
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")

  #Get predicted values
  psi <- predict(object, 'psi', na.rm=FALSE)$Predicted
  lam <- predict(object, 'det', na.rm=FALSE)$Predicted
  if(T>1){
    p_col <- predict(object, 'col', na.rm=FALSE)$Predicted
    p_col <- matrix(p_col, N, T, byrow=TRUE)
    p_ext <- predict(object, 'ext', na.rm=FALSE)$Predicted
    p_ext <- matrix(p_ext, N, T, byrow=TRUE)
  }

  tmax <- object@data@surveyLength
  not_na <- which(!is.na(lam))

  simlist <- list()
  for(s in 1:nsim){
    ttd <- rep(NA, length(lam))
    if(tdist == "weibull"){
      k <- exp(coef(object)['k(k)'])
      ttd[not_na] <- stats::rweibull(length(not_na),k,1/lam[not_na])
    } else {
      ttd[not_na] <- stats::rexp(length(not_na), lam[not_na])
    }
    #Truncate
    ttd <- matrix(ttd, nrow=N, byrow=T)
    ttd[which(ttd>tmax)] <- tmax[which(ttd>tmax)]

    #Latent state
    z <- matrix(NA, N, T)
    z[,1] <- rbinom(N, 1, psi)

    if(T>1){
      for (t in 1:(T-1)){
        z_ext <- rbinom(N, 1, 1-p_ext[,t])
        z_col <- rbinom(N, 1, p_col[,t])
        z[,t+1] <- ifelse(z[,t], z_ext, z_col)
      }
    }

    #Detection process
    yout <- matrix(NA, N, J*T)
    d_ind <- 1
    for (t in 1:T){
      for (j in 1:J){
        yout[,d_ind] <- ifelse(z[,t], ttd[,d_ind], tmax[,d_ind])
        d_ind <- d_ind + 1
      }
    }

    #Add NAs
    nas <- which(is.na(object@data@y))
    yout[nas] <- NA
    simlist[[s]] <- yout
  }
  simlist
})

setMethod("simulate", "unmarkedFitNmixTTD",
          function(object,  nsim = 1, seed = NULL, na.rm = FALSE)
{

  M <- nrow(object@data@y)
  J <- ncol(object@data@y)
  tdist <- ifelse("shape" %in% names(object@estimates), "weibull", "exp")
  mix <- ifelse("alpha" %in% names(object@estimates), "NB", "P")

  #Get predicted values
  abun <- predict(object, 'state', na.rm=FALSE)$Predicted
  lam <- predict(object, 'det', na.rm=FALSE)$Predicted
  tmax <- as.vector(t(object@data@surveyLength))
  not_na <- !is.na(lam)

  simlist <- list()
  for(s in 1:nsim){

    if(mix=="P"){
      N <- rpois(M, abun)
    } else if(mix=="NB"){
      alpha <- exp(coef(object, "alpha"))
      N <- rnbinom(M, mu=abun, size=alpha)
    }

    lamN <- lam*rep(N, each=J)
    is_zero <- lamN==0
    to_sim <- !is_zero & not_na

    ttd <- rep(NA, length(lamN))
    ttd[is_zero] <- tmax[is_zero]

    if(tdist == "weibull"){
      k <- exp(coef(object)['k(k)'])
      ttd[to_sim] <- stats::rweibull(sum(to_sim),k,1/lamN[to_sim])
    } else {
      ttd[to_sim] <- stats::rexp(sum(to_sim), lamN[to_sim])
    }
    #Truncate
    ttd[ttd>tmax] <- tmax[ttd>tmax]

    ttd <- matrix(ttd, nrow=M, byrow=T)

    #Add NAs
    nas <- which(is.na(object@data@y))
    ttd[nas] <- NA
    simlist[[s]] <- ttd
  }
  simlist
})

setMethod("simulate", "unmarkedFitColExt",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    data <- object@data
    psiParms <- coef(object, 'psi')
    detParms <- coef(object, 'det')
    colParms <- coef(object, 'col')
    extParms <- coef(object, 'ext')
    formulaList <- list(psiformula=object@psiformula,
        gammaformula=object@gamformula,
        epsilonformula=object@epsformula,
        pformula=object@detformula)
    designMats <- getDesign(object@data, object@formula)
    V.itj <- designMats$V
    X.it.gam <- designMats$X.gam
    X.it.eps <- designMats$X.eps
    W.i <- designMats$W
    y <- designMats$y

    M <- nrow(y)	# M <- nrow(X.it)
    nY <- data@numPrimary
    J <- obsNum(data)/nY

    psiP <- plogis(W.i %*% psiParms)
    detP <- plogis(V.itj %*% detParms)
    colP <- plogis(X.it.gam  %*% colParms)
    extP <- plogis(X.it.eps %*% extParms)

    detP <- array(detP, c(J, nY, M))
    detP <- aperm(detP, c(3, 1, 2))
    colP <- matrix(colP, M, nY, byrow = TRUE)
    extP <- matrix(extP, M, nY, byrow = TRUE)

    simList <- list()
    for(s in 1:nsim) {
        ## generate first year's data
        x <- matrix(0, M, nY)
        x[,1] <- rbinom(M, 1, psiP)

        ## create transition matrices (phi^T)
        phis <- array(NA,c(2,2,nY-1,M)) #array of phis for each
        for(i in 1:M) {
            for(t in 1:(nY-1)) {
                phis[,,t,i] <- matrix(c(1-colP[i,t], colP[i,t], extP[i,t],
                    1-extP[i,t]))
                }
            }

        ## generate latent years 2:T
        for(i in 1:M) {
            for(t in 2:nY) {
                x[i,t] <- rbinom(1, 1, phis[2,x[i,t-1]+1,t-1,i])
                }
            }

        ## generate observations
        y <- array(NA, c(M, J, nY))
        for(t in 1:nY) {
            y[,,t] <- rbinom(M*J, 1, x[,t]*detP[,,t])
            }

        y.mat <- y[,,1]
        for(i in 2:dim(y)[3]) {
            y.mat <- cbind(y.mat,y[,,i])
            }
        simList[[s]] <- y.mat
        }
    return(simList)
})




setMethod("simulate", "unmarkedFitOccuRN",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    formula <- object@formula
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y; X <- designMats$X; V <- designMats$V
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    detParms <- coef(object, 'det')
    r.ij <- plogis(V %*% detParms)
    r <- matrix(r.ij, M, J, byrow = TRUE)
    lamParms <- coef(object, 'state')
    lambda <- exp(X %*% lamParms + X.offset)
    simList <- list()
    for(s in 1:nsim) {
        N.i <- rpois(M, lambda)
        N.ij <- rep(N.i, each = J)
        y <- matrix(NA, M, J)
        for(i in 1:J) {
            y[,i] <- rbinom(M, N.i, r[,i])
            }
        simList[[s]] <- ifelse(y > 0, 1, 0)
        }
    return(simList)
})



setMethod("simulate", "unmarkedFitGMM",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    formula <- object@formula
    umf <- object@data
    mixture <- object@mixture
    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xlam <- D$Xlam
    Xphi <- D$Xphi
    Xdet <- D$Xdet

    Xlam.offset <- D$Xlam.offset
    Xphi.offset <- D$Xphi.offset
    Xdet.offset <- D$Xdet.offset
    if (is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
    if (is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
    if (is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

    n <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T

    lamParms <- coef(object, type = "lambda")
    phiParms <- coef(object, type = "phi")
    detParms <- coef(object, type = "det")
    lam <- drop(exp(Xlam %*% lamParms + Xlam.offset))
    if(is.null(phiParms))
        phi <- rep(1, nrow(Xphi))
    else
        phi <- as.numeric(plogis(Xphi %*% phiParms + Xphi.offset))
    phi.mat <- matrix(phi, nrow=n, ncol=T, byrow=TRUE)
    p <- as.numeric(plogis(Xdet %*% detParms + Xdet.offset))

    cp.arr <- array(NA, c(n, T, J+1))
    cp.mat <- getP(object, na.rm = na.rm)
    cp.mat[is.na(y)] <- NA
    cp.temp <- array(cp.mat, c(n, J, T))
    cp.arr[,,1:J] <- aperm(cp.temp, c(1,3,2))
    cp.arr[,, 1:J][is.na(y)]<- NA   # Andy added 5/30
    cp.arr[,,J+1] <- 1 - apply(cp.arr[,,1:J,drop=FALSE], 1:2, sum, na.rm=TRUE)

    simList <- list()
    for(s in 1:nsim) {
        switch(mixture,
            P = M <- rpois(n=n, lambda=lam),
            NB = M <- rnbinom(n=n, mu=lam,
                size=exp(coef(object, type="alpha"))),
            ZIP = {
              psi <- plogis(coef(object['psi']))
              M <- rzip(n, lambda=lam, psi=psi)
            }
        )

        N <- rbinom(n*T, size=M, prob=phi.mat)
        # bug fix 3/16/2010
        N <- matrix(N, nrow=n, ncol=T, byrow=FALSE) # , byrow=TRUE)

        y.sim <- array(NA, c(n, J, T))
        for(i in 1:n) {
            for(t in 1:T) {
                if(is.na(N[i,t]))
                    next
                pi.it <- cp.arr[i,t,]
                na.it <- is.na(pi.it)
                pi.it[na.it] <- 0
                y.sim[i,,t] <- drop(rmultinom2(1, N[i,t], pi.it))[1:J]
                y.sim[i,na.it[1:J],t] <- NA
            }
        }
        simList[[s]] <- matrix(y.sim, nrow=n, ncol=J*T) # note, byrow=F
        }
    return(simList)
})





setMethod("simulate", "unmarkedFitGPC",
    function(object, nsim = 1, seed = NULL, na.rm = TRUE)
{
    formula <- object@formula
    umf <- object@data
    mixture <- object@mixture
    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xlam <- D$Xlam
    Xphi <- D$Xphi
    Xdet <- D$Xdet

    Xlam.offset <- D$Xlam.offset
    Xphi.offset <- D$Xphi.offset
    Xdet.offset <- D$Xdet.offset
    if (is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
    if (is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
    if (is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

    R <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T

    lamParms <- coef(object, type = "lambda")
    phiParms <- coef(object, type = "phi")
    detParms <- coef(object, type = "det")
    lam <- drop(exp(Xlam %*% lamParms + Xlam.offset))
    if(is.null(phiParms))
        phi <- rep(1, nrow(Xphi))
    else
        phi <- as.numeric(plogis(Xphi %*% phiParms + Xphi.offset))
    phi.mat <- matrix(phi, nrow=R, ncol=T, byrow=TRUE)
    p <- as.numeric(plogis(Xdet %*% detParms + Xdet.offset))
    p <- matrix(p, R, J*T, byrow=TRUE)
    p <- array(p, c(R, J, T))

    simList <- list()
    for(s in 1:nsim) {
        switch(mixture,
            P = M <- rpois(n=R, lambda=lam),
               #               FIXME: Add ZIP
            NB = M <- rnbinom(n=R, mu=lam,
               size=exp(coef(object, type="alpha"))),   
            ZIP = {
              psi <- plogis(coef(object['psi']))
              M <- rzip(R, lambda=lam, psi=psi)
            }
        )

        N <- rbinom(R*T, size=M, prob=phi.mat)
        N <- matrix(N, nrow=R, ncol=T, byrow=FALSE)

        y.sim <- array(NA, c(R, J, T))
        for(i in 1:R) {
            for(t in 1:T) {
                if(is.na(N[i,t]))
                    next
                y.sim[i,,t] <- rbinom(J, N[i,t], p[i,,t])
            }
        }
        y.sim <- matrix(y.sim, nrow=R, ncol=J*T)
        y.sim[is.na(y)] <- NA # Not necessary if covariates exist!
        simList[[s]] <- y.sim
    }
    return(simList)
})


setMethod("simulate", "unmarkedFitGDS",
    function(object, nsim = 1, seed = NULL, na.rm=TRUE)
{
    formula <- object@formula
    umf <- object@data
    db <- umf@dist.breaks
    w <- diff(db)
    mixture <- object@mixture
    keyfun <- object@keyfun

    D <- getDesign(umf, formula, na.rm = na.rm)
    y <- D$y
    Xlam <- D$Xlam
    Xphi <- D$Xphi
    Xdet <- D$Xdet
    Xlam.offset <- D$Xlam.offset
    Xphi.offset <- D$Xphi.offset
    Xdet.offset <- D$Xdet.offset

    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
    if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
    if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

    M <- nrow(y)
    T <- umf@numPrimary
    R <- ncol(y)
    J <- R / T

    lamPars <- coef(object, type="lambda")
    if(T>1)
        phiPars <- coef(object, type="phi")
    detPars <- coef(object, type="det")

    lambda <- exp(Xlam %*% lamPars + Xlam.offset)
    if(T==1)
        phi <- matrix(1, M, T)
    else {
        phi <- plogis(Xphi %*% phiPars + Xphi.offset)
        phi <- matrix(phi, M, T, byrow=TRUE)
        }

    if(identical(object@output, "density")) {
        switch(umf@survey,
            line = {
                tlength <- umf@tlength
                A <- tlength * max(db) * 2
                },
            point = {
                A <- pi * max(db)^2
                })
        switch(umf@unitsIn,
            m = A <- A / 1e6,
            km = A <- A)
        switch(object@unitsOut,
            ha = A <- A * 100,
            kmsq = A <- A)
        lambda <- lambda * A
        }
    cp <- getP(object, na.rm = na.rm)
    ysim <- cpa <- array(cp, c(M, J, T))

    simList <- list()
    for(s in 1:nsim) {
        for(i in 1:M) {
            switch(mixture,
                P = Ns <- rpois(1, lambda[i]),
                NB = {
                    alpha <- exp(coef(object, type="alpha"))
                    Ns <- rnbinom(1, mu=lambda[i], size=alpha)
                    },
                ZIP = {
                    psi <- plogis(coef(object['psi']))
                    Ns <- rzip(1, lambda[i], psi)
                    }
            )
            for(t in 1:T) {
                N <- rbinom(1, Ns, phi[i,t])
                cp.it <- cpa[i,,t]
                cp.it[J+1] <- 1-sum(cp.it)
                y.it <- as.integer(rmultinom2(1, N, prob=cp.it))
                ysim[i,,t] <- y.it[1:J]
                }
            }
        y.mat <- matrix(ysim, nrow=M)
        simList[[s]] <- y.mat
        }
    return(simList)
})



