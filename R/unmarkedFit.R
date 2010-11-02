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
        covMatBS = "optionalMatrix")) # list of bootstrap sample fits

# constructor for unmarkedFit objects
unmarkedFit <- function(fitType, call, formula, data, sitesRemoved, estimates, 
    AIC, opt, negLogLike, nllFun) 
{
    umfit <- new("unmarkedFit", fitType = fitType, call = call, 
        formula = formula, data = data, sitesRemoved = sitesRemoved, 
        estimates = estimates, AIC = AIC, opt = opt, negLogLike = negLogLike, 
        nllFun = nllFun)
    return(umfit)
}

# ---------------------------- CHILD CLASSES -----------------------------------


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


setClass("unmarkedFitOccu", 
    representation(knownOcc = "logical"),
    contains = "unmarkedFit")			


setClass("unmarkedFitMPois", 
    contains = "unmarkedFit")			


setClass("unmarkedFitOccuRN", 
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
        

# -------------------------- Show and Summary ----------------------------------


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
        cat("\nSites removed:", object@sitesRemoved)
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
    function(obj, coefficients, type, offset = NULL) 
{
    stopifnot(!missing(type))
    stopifnot(type %in% names(obj))
    estimate <- obj@estimates[type]
    linearComb(estimate, coefficients, offset)
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

setMethod("[", "unmarkedFit", function(x, i, j, drop) 
{
    x@estimates[i]
})

setMethod("names", "unmarkedFit", function(x) 
{
    names(x@estimates)
})



# ----------------------------- Prediction ------------------------------------



setMethod("predict", "unmarkedFit", 
    function(object, type, newdata, backTransform = TRUE, na.rm = TRUE, 
        appendData = FALSE, ...) 
{
    if(missing(newdata) || is.null(newdata))
        newdata <- getData(object)
    formula <- object@formula
    detformula <- as.formula(formula[[2]])
    stateformula <- as.formula(paste("~", formula[3], sep=""))
    if(inherits(newdata, "unmarkedFrame"))
        class(newdata) <- "unmarkedFrame"
    cls <- class(newdata)
    switch(cls, 
    unmarkedFrame = {
        designMats <- getDesign(newdata, formula, na.rm = na.rm)
        switch(type, 
            state = {
                X <- designMats$X
                offset <- designMats$X.offset
                },
            det = {
                X <- designMats$V
                offset <- designMats$V.offset
                })
        },
    data.frame = {
        switch(type, 
            state = {
                mf <- model.frame(stateformula, newdata)
                X <- model.matrix(stateformula, mf)
                offset <- model.offset(mf)
                },
            det = {
                mf <- model.frame(detformula, newdata)
                X <- model.matrix(detformula, mf)
                offset <- model.offset(mf)
                })
            })
    out <- data.frame(matrix(NA, nrow(X), 2,  
        dimnames=list(NULL, c("Predicted", "SE"))))
    lc <- linearComb(object, X, type, offset = offset)
    if(backTransform) lc <- backTransform(lc)
    out$Predicted <- coef(lc)
    out$SE <- SE(lc)
    ci <- as.data.frame(confint(lc))
    colnames(ci) <- c("lower", "upper")
    out <- cbind(out, ci)
    if(appendData)
        out <- data.frame(out, as(newdata, "data.frame"))
    return(out)
})



setMethod("predict", "unmarkedFitColExt", 
    function(object, type, newdata, backTransform = TRUE, na.rm = TRUE, 
        appendData = FALSE, ...) 
{
    if(missing(newdata) || is.null(newdata))
        newdata <- getData(object)
    formula <- object@formula
    cls <- class(newdata)
    switch(cls, 
    unmarkedMultFrame = {
        designMats <- getDesign(newdata, formula, na.rm = na.rm)
        switch(type, 
            psi = {
                X <- designMats$W
                #offset <- designMats$W.offset
                },
            col = X <- designMats$X.gam,
            ext = X <- designMats$X.eps,                
            det = {
                X <- designMats$V
                #offset <- designMats$V.offset
                })
            },
    data.frame = {
        aschar1 <- as.character(formula)
        aschar2 <- as.character(formula[[2]])
        aschar3 <- as.character(formula[[2]][[2]])    

        detformula <- as.formula(paste(aschar1[1], aschar1[3]))
        epsformula <- as.formula(paste(aschar2[1], aschar2[3]))
        gamformula <- as.formula(paste(aschar3[1], aschar3[3]))
        psiformula <- as.formula(formula[[2]][[2]][[2]])
            
        switch(type, 
            psi = {
                mf <- model.frame(psiformula, newdata)
                X <- model.matrix(psiformula, mf)
                #offset <- model.offset(mf)
                },
            col = {
                mf <- model.frame(gamformula, newdata)
                X <- model.matrix(gamformula, mf)
                #offset <- model.offset(mf)
                },
            ext = {
                mf <- model.frame(epsformula, newdata)
                X <- model.matrix(epsformula, mf)
                #offset <- model.offset(mf)
                },               
            det = {
                mf <- model.frame(detformula, newdata)
                X <- model.matrix(detformula, mf)
                #offset <- model.offset(mf)
                })
            })
    out <- data.frame(matrix(NA, nrow(X), 2, 
        dimnames=list(NULL, c("Predicted", "SE"))))
    lc <- linearComb(object, X, type)#, offset = offset)
    if(backTransform) lc <- backTransform(lc)
    out$Predicted <- coef(lc)
    out$SE <- SE(lc)
    ci <- as.data.frame(confint(lc))
    colnames(ci) <- c("lower", "upper")
    out <- cbind(out, ci)
    if(appendData)
        out <- data.frame(out, as(newdata, "data.frame"))
    return(out)
})









setMethod("predict", "unmarkedFitGMM", 
    function(object, type, newdata, backTransform = TRUE, na.rm = TRUE, 
        appendData = FALSE, ...) 
{
    if(missing(newdata) || is.null(newdata))
        newdata <- getData(object)
    formlist <- object@formlist
    lamformula <- formlist$lambdaformula
    phiformula <- formlist$phiformula        
    pformula <- formlist$pformula        
    formula <- object@formula
        
    if(inherits(newdata, "unmarkedFrame"))
        cls <- "unmarkedFrame"
    switch(cls, 
        unmarkedFrame = {
            D <- getDesign(newdata, formula, na.rm = na.rm)
            switch(type, 
                lambda = {
                    X <- D$Xlam
                    offset <- D$Xlam.offset
                    },
                phi = {
                    X <- D$Xphi
                    offset <- D$Xphi.offset
                    },
                det = {   # Note, this is p not pi
                    X <- D$Xdet
                    offset <- D$Xdet.offset
                    })
              },
        data.frame = {
            switch(type, 
                lambda = {
                    mf <- model.frame(lambdaformula, newdata)
                    X <- model.matrix(lambdaformula, mf)
                    offset <- model.offset(mf)
                    },
                phi = {
                    mf <- model.frame(phiformula, newdata)
                    X <- model.matrix(phiformula, mf)
                    offset <- model.offset(mf)
                    },                
                det = {   # Note, this is p not pi
                  mf <- model.frame(detformula, newdata)
                  X <- model.matrix(detformula, mf)
                  offset <- model.offset(mf)
                })
            })
        out <- data.frame(matrix(NA, nrow(X), 2, 
            dimnames=list(NULL, c("Predicted", "SE"))))
        lc <- linearComb(object, X, type, offset = offset)
        if(backTransform) lc <- backTransform(lc)
        out$Predicted <- coef(lc)
        out$SE <- SE(lc)
        ci <- as.data.frame(confint(lc))
        colnames(ci) <- c("lower", "upper")
        out <- cbind(out, ci)
        if(appendData)
            out <- data.frame(out, as(newdata, "data.frame"))
        return(out)
        })




# ---------------------- coef, vcov, and SE -----------------------------------


setMethod("coef", "unmarkedFit",
    function(object, type, altNames = TRUE) 
{
    if(missing(type)) {
        co <- lapply(object@estimates@estimates, 
            function(x) coef(x, altNames=altNames))
        names(co) <- NULL
        co <- unlist(co)
    } else {
        co <- coef(object[type], altNames=altNames)
        }
    co
})


setMethod("vcov", "unmarkedFit",
    function (object, type, altNames = TRUE, method = "hessian", ...) 
{
    method <- match.arg(method, c("hessian", "nonparboot"))
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
        })
    rownames(v) <- colnames(v) <- names(coef(object, altNames=altNames))
    if (missing(type)) {
        return (v)
    } else {
        inds <- .estimateInds(object)[[type]]
        return (v[inds, inds, drop = FALSE])
        }
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
    chisq <- 2 * diff(c(logLik(m1), logLik(m2)))
    DF <- diff(c(length(coef(m1)), length(coef(m2))))
    pval <- pchisq(chisq, DF, lower=FALSE)
    return(data.frame(Chisq=chisq, DF = DF, 'Pr(>Chisq)' = pval, check.names=F))
}) 
    
    



setMethod("confint", "unmarkedFit", function(object, parm, level = 0.95, type, 
    method = c("normal", "profile")) 
{
    method <- match.arg(method)
    if(missing(type)) 
        stop(paste("Must specify type as one of (", paste(names(object),collapse=", "),").",sep=""))
    if(missing(parm)) 
        parm <- 1:length(object[type]@estimates)
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
    X <- des$X
    X.offset <- des$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    state <- do.call(object['state']@invlink, 
        list(X %*% coef(object, 'state') + X.offset))
    state <- as.numeric(state)  ## E(X) for most models
    p <- getP(object, na.rm = na.rm) # P(detection | presence)
    fitted <- state * p  # true for models with E[Y] = p * E[X]
    fitted
})		



setMethod("fitted", "unmarkedFitDS", function(object, na.rm = FALSE) 
{
    data <- object@data
    D <- getDesign(data, object@formula, na.rm = na.rm)
    X <- D$X
    X.offset <- D$X.offset
    if (is.null(X.offset)) {
      X.offset <- rep(0, nrow(X))
    }
    lambda <- drop(exp(X %*% coef(object, 'state') + X.offset))
    a <- calcAreas(dist.breaks = data@dist.breaks, tlength = data@tlength, 
	   survey = data@survey, output = object@output, M = numSites(data), 
	   J = ncol(getY(data)), unitsIn = data@unitsIn, unitsOut = object@unitsOut)
    if(length(D$removed.sites)>0)
        a <- a[-D$removed.sites,]
    p <- getP(object, na.rm = na.rm)
    fitted <- lambda * p * a
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
    X <- des$X
    X.offset <- des$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    y <- des$y	# getY(data) ... to be consistent w/NA handling?
    M <- nrow(X)
    J <- ncol(y)
    state <- exp(X %*% coef(object, 'state') + X.offset)
    p <- getP(object, na.rm = na.rm)
    mix <- object@mixture
    switch(mix,
        P = {
            fitted <- as.numeric(state) * p 
            },
        NB = {
            if(missing(K)) K <- max(y, na.rm = TRUE) + 20 
            k <- 0:K
            k.ijk <- rep(k, M*J)
            state.ijk <- state[rep(1:M, each = J*(K+1))]
            alpha <- exp(coef(object['alpha']))
            prob.ijk <- dnbinom(k.ijk, mu = state.ijk, size = alpha)
            all <- cbind(rep(as.vector(t(p)), each = K + 1), k.ijk, prob.ijk)
            prod.ijk <- rowProds(all)
            fitted <- colSums(matrix(prod.ijk, K + 1, M*J))
            fitted <- matrix(fitted, M, J, byrow = TRUE)
            })
    return(fitted)
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
    if(missing(K)) K <- max(y, na.rm = TRUE) + 20 
            
    lam <- rep(lam, each = J)
            
    fitted <- 1 - exp(-lam*r) ## analytical integration.
            
    return(matrix(fitted, M, J, byrow = TRUE))
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
    designMats <- getDesign(object@data, formlist = formulaList)
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
    phi <- plogis(Xphi %*% coef(object, 'phi') + Xphi.offset)
    phi.mat <- matrix(phi, nrow=M, ncol=T, byrow=TRUE)
    phi.ijt <- as.numeric(apply(phi, 2, rep, times=J))
    cp <- getP(object, na.rm = na.rm)

# My tests show that E[y] is the same under the P and NB models    
#    mix <- object@mixture
#    switch(mix,
#        P = {
#            fitted <- lambda * phi.ijt * as.numeric(cp) # recycle
#            fitted <- matrix(fitted, M, J*T) 
#            },
#        NB = {
#            K <- object@K 
#            k <- 0:K
#            k.ijk <- rep(k, M*J)
#            lambda.ijk <- lambda[rep(1:M, each = J*(K+1))]
#            alpha <- exp(coef(object['alpha']))
#            prob.ijk <- dnbinom(k.ijk, mu = lambda.ijk, size = alpha)
#            all <- cbind(rep(as.vector(t(cp)), each = K + 1), k.ijk, prob.ijk)
#            prod.ijk <- rowProds(all)
#            fitted <- colSums(matrix(prod.ijk, K + 1, M*J))
#            fitted <- matrix(fitted, M, J*T, byrow = TRUE)
#            })

    fitted <- lambda * phi.ijt * as.numeric(cp) # recycle
    fitted <- matrix(fitted, M, J*T)
    return(fitted)
})
          



setMethod("profile", "unmarkedFit",
          function(fitted, type, parm, seq) {
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
          function(object) {
            object@opt$hessian
          })


setMethod("update", "unmarkedFit", 
    function(object, formula., ..., evaluate = TRUE) 
{
    call <- object@call
    origFormula <- formula(call)
    if (is.null(call)) 
        stop("need an object with call slot")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(formula.)) {
        detformula <- as.formula(origFormula[[2]])
        stateformula <- as.formula(paste("~", origFormula[3], sep=""))
        newDetformula <- as.formula(formula.[[2]])
        upDetformula <- update.formula(detformula, newDetformula)
        newStateformula <- as.formula(paste("~", formula.[3], sep=""))
        upStateformula <- update.formula(stateformula, newStateformula)
        call$formula <- as.formula(paste(
			deparse(upDetformula, width=500), 
            deparse(upStateformula, width=500)))
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
        eval(call, parent.frame())
    else call
})


setMethod("update", "unmarkedFitColExt", 
    function(object, ..., evaluate = TRUE) 
{
    call <- object@call
    if (is.null(call)) 
        stop("need an object with call slot")
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
            if (any(!existing)) {
                call <- c(as.list(call), extras[!existing])
                call <- as.call(call)
            }
        }
    if (evaluate) 
        eval(call, parent.frame())
    else call
})
          
          

setMethod("update", "unmarkedFitGMM", 
    function(object, lambdaformula, phiformula, pformula, ..., evaluate = TRUE)
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
    extras <- match.call(expand.dots = FALSE)$...
    if(length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
            }
        }
    if (evaluate) 
        eval(call, parent.frame())
    else call
})

          
          

setGeneric("sampleSize", function(object) standardGeneric("sampleSize"))
setMethod("sampleSize", "unmarkedFit",
          function(object) {
            data <- getData(object)
            M <- numSites(data)
            M <- M - length(object@sitesRemoved)
            M
          })


setGeneric("getData", function(object) standardGeneric("getData"))	
setMethod("getData", "unmarkedFit",
          function(object) {
            object@data
          })


setGeneric("nllFun", function(object) standardGeneric("nllFun"))
setMethod("nllFun", "unmarkedFit", function(object) object@nllFun)

setGeneric("mle", function(object) standardGeneric("mle"))
setMethod("mle", "unmarkedFit", function(object) object@opt$par)

setClass("profile", representation(prof = "matrix"))

setGeneric("smoothed", function(object, mean=TRUE) standardGeneric("smoothed"))
setMethod("smoothed","unmarkedFitColExt",
          function(object, mean) {
            if(mean) object@smoothed.mean
            else object@smoothed
          })

setGeneric("projected", function(object, mean=TRUE) standardGeneric("projected"))
setMethod("projected","unmarkedFitColExt",
          function(object, mean) {
            if(mean) object@projected.mean
            else object@projected
          })

setMethod("plot", c("profile", "missing"),
          function(x) {
            plot(x@prof[,1], x@prof[,2], type = "l")
          })


setMethod("residuals", "unmarkedFit", function(object, ...)
          {
            y <- getY(object@data)
            e <- fitted(object, na.rm = FALSE)	 
            r <- y - e
            return(r)
          })

setMethod("residuals", "unmarkedFitOccu", function(object, ...)
          {
            y <- getY(object@data)
            y <- truncateToBinary(y)
            e <- fitted(object, na.rm = FALSE)	 
            r <- y - e
            return(r)
          })

setMethod("residuals", "unmarkedFitOccuRN", function(object, ...)
          {
            y <- getY(object@data)
            y <- truncateToBinary(y)
            e <- fitted(object, na.rm = FALSE)	 
            r <- y - e
            return(r)
          })


setMethod("plot", c(x = "unmarkedFit", y = "missing"), 
          function(x, y, ...)
          {
            r <- residuals(x)
            e <- fitted(x, na.rm = FALSE)
            plot(e, r, ylab = "Residuals", xlab = "Predicted values")
            abline(h = 0, lty = 3, col = "gray")
          })




setMethod("hist", "unmarkedFitDS", 
    function(x, lwd=1, lty=1, ...) {
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
                int <- 2 * integrate(dnorm, dbreaks[1], dbreaks[nb], 
                    sd=sigma)$value
                h$density <- h$density * int
                plot(h, freq=F, ...)
                plot(function(x) 2 * dnorm(x, mean=0, sd=sigma), 
                    min(dbreaks), max(dbreaks), add=T, lwd=lwd, lty=lty)
                },
            point = {
                int <- integrate(drhn, dbreaks[1], dbreaks[nb], 
                    sigma=sigma)$value
                h$density <- h$density * int
                plot(h, freq=F, ...)
                plot(function(r) drhn(r, sigma=sigma), min(dbreaks), 
                    max(dbreaks), add=T, lwd=lwd, lty=lty)
                })
            },
        exp = {		# This doesn't work on example fm4
            rate <- exp(coef(x, type="det"))
            if(length(rate) > 1)
                stop("This method only works when there are no detection covars")
            switch(survey,
            line = {
                int <- integrate(dxexp, dbreaks[1], dbreaks[nb], rate=rate)$value
                h$density <- h$density * int
                plot(h, freq=F, ...)
                plot(function(x) dxexp(x, rate=rate), min(dbreaks), 
                    max(dbreaks), add=T, lwd=lwd, lty=lty)
                },
            point = {
                int <- integrate(drexp, dbreaks[1], dbreaks[nb], rate=rate)$value
                h$density <- h$density * int
                plot(h, freq=F, ...)
                plot(function(r) drexp(r, rate=rate), min(dbreaks), 
                    max(dbreaks), add=T, lwd=lwd, lty=lty)
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
                plot(function(x) dxhaz(x, shape=shape, scale=scale), 
                    min(dbreaks), max(dbreaks), add=T, lwd=lwd, lty=lty)
                },
            point = {
                int <- integrate(drhaz, dbreaks[1], dbreaks[nb], 
                    shape=shape, scale=scale)$value
                    h$density <- h$density * int
                    plot(h, freq=F, ...)
                    plot(function(r) drhaz(r, shape=shape, scale=scale), 
                        min(dbreaks), max(dbreaks), add=T, lwd=lwd, lty=lty)
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
                    min(dbreaks), max(dbreaks), add=T, lwd=lwd, lty=lty)
                }
            )}
            )
        })




############################# CHILD CLASS METHODS ##############################

                                        # Extract detection probs
setGeneric("getP", function(object, ...) standardGeneric("getP"))


setMethod("getP", "unmarkedFit", function(object, na.rm = TRUE) 
          {
            formula <- object@formula
            detformula <- as.formula(formula[[2]])
            umf <- object@data
            designMats <- getDesign(umf, formula, na.rm = na.rm)
            y <- designMats$y
            V <- designMats$V
            V.offset <- designMats$V.offset
            if (is.null(V.offset)) {
              V.offset <- rep(0, nrow(V))
            }
            M <- nrow(y)
            J <- ncol(y)
            ppars <- coef(object, type = "det")
            p <- plogis(V %*% ppars + V.offset)
            p <- matrix(p, M, J, byrow = TRUE)
            return(p)
          })




setMethod("getP", "unmarkedFitDS", 
    function(object, na.rm = TRUE) {
        formula <- object@formula
        detformula <- as.formula(formula[[2]])
        umf <- object@data
        designMats <- getDesign(umf, formula, na.rm = na.rm)
        y <- designMats$y
        V <- designMats$V
        V.offset <- designMats$V.offset
        if (is.null(V.offset)) {
          V.offset <- rep(0, nrow(V))
        }
        M <- nrow(y)
        J <- ncol(y)
        ppars <- coef(object, type = "det")
        d <- umf@dist.breaks
        survey <- umf@survey
        key <- object@keyfun
        switch(key, 
        halfnorm = {
            sigma <- exp(V %*% ppars + V.offset)
            p <- sapply(sigma, function(x) cp.hn(d = d, s = x, survey = survey))
            }, 
        exp = {
            rate <- exp(V %*% ppars + V.offset)
            p <- sapply(rate, function(x) cp.exp(d = d, r = x, survey = survey))
            }, 
        hazard = {
            shape <- exp(V %*% ppars + V.offset)
            scale <- exp(coef(object, type="scale"))
            p <- sapply(shape, function(x) cp.haz(d = d, shape = x, 
            scale = scale, survey = survey))
            },
		uniform = p <-1)
        p <- matrix(p, M, J, byrow = TRUE)
        return(p)
        })






setMethod("getP", "unmarkedFitMPois", function(object, na.rm = TRUE) 
          {
            formula <- object@formula
            detformula <- as.formula(formula[[2]])
            piFun <- object@data@piFun
            umf <- object@data
            designMats <- getDesign(umf, formula, na.rm = na.rm)
            y <- designMats$y
            V <- designMats$V
            V.offset <- designMats$V.offset
            if (is.null(V.offset)) {
              V.offset <- rep(0, nrow(V))
            }
            M <- nrow(y)
            J <- ncol(y)
            ppars <- coef(object, type = "det")
            p <- plogis(V %*% ppars + V.offset)
            p <- matrix(p, M, J, byrow = TRUE)
            pi <- do.call(piFun, list(p = p))
            return(pi)
          })


setMethod("getP", "unmarkedFitColExt", function(object, na.rm = TRUE)
          {
            stop("getP is not yet implemented for colext fits.")
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
    R <- ncol(y) 
    J <- R / T
    
    ppars <- coef(object, type = "det")
    p <- plogis(Xdet %*% ppars + Xdet.offset)
    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, J, T))
    p <- aperm(p, c(1,3,2))     

    cp <- array(as.numeric(NA), c(M, T, J))
    for(t in 1:T) cp[,t,] <- do.call(piFun, list(p[,t,]))
    cp <- aperm(cp, c(1,3,2))
    cp <- matrix(cp, nrow=M, ncol=R)
    
    return(cp)
})





setMethod("simulate", "unmarkedFitDS", 
    function(object, nsim = 1, seed = NULL, na.rm=TRUE)
    {
    formula <- object@formula
    umf <- object@data
    designMats <- getDesign(umf, formula, na.rm = na.rm)
    y <- designMats$y
    X <- designMats$X
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
      X.offset <- rep(0, nrow(X))
    }
    a <- calcAreas(dist.breaks = umf@dist.breaks, tlength = umf@tlength, 
	   survey = umf@survey, output = object@output, M = numSites(umf), 
	   J = ncol(getY(umf)), unitsIn = umf@unitsIn, unitsOut = object@unitsOut)
    if(length(designMats$removed.sites)>0)
        a <- a[-designMats$removed.sites,]
    M <- nrow(y)
    J <- ncol(y)
    lamParms <- coef(object, type = "state")
    lam <- drop(exp(X %*% lamParms + X.offset))
    pmat <- getP(object, na.rm = na.rm)
    simList <- list()
    for(i in 1:nsim) {
        yvec <- rpois(M * J, lam * pmat * a)
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
    X <- designMats$X
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    allParms <- coef(object, altNames = FALSE)
    lamParms <- coef(object, type = "state")
    lam <- as.numeric(exp(X %*% lamParms + X.offset)) 
    lamvec <- rep(lam, each = J)
    pvec <- c(t(getP(object, na.rm = na.rm)))
    mix <- object@mixture
    simList <- list()
    for(i in 1:nsim) {
        switch(mix, 
            P = yvec <- rpois(M * J, lamvec * pvec),
            NB = {
                N <- rnbinom(M, size = exp(coef(object["alpha"])), mu = lam)
                    yvec <- rbinom(M * J, size = rep(N, each = J), prob = pvec)
                })
        simList[[i]] <- matrix(yvec, M, J, byrow = TRUE)
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
    X <- designMats$X
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    lamParms <- coef(object, type = "state")
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
    X <- designMats$X
    X.offset <- designMats$X.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
        }
    M <- nrow(y)
    J <- ncol(y)
    allParms <- coef(object, altNames = FALSE)
    psiParms <- coef(object, type = "state")
    psi <- as.numeric(plogis(X %*% psiParms + X.offset))
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
    designMats <- getDesign(object@data, formlist = formulaList)
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
    designMats <- unmarked:::getDesign(umf, formula, na.rm = na.rm)
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
    phi <- as.numeric(plogis(Xphi %*% phiParms + Xphi.offset))
    phi.mat <- matrix(phi, nrow=n, byrow=TRUE)
    p <- as.numeric(plogis(Xdet %*% detParms + Xdet.offset))    

    cp.arr <- array(NA, c(n, T, J+1))  
    cp.mat <- getP(object, na.rm = na.rm)
    cp.temp <- array(cp.mat, c(n, J, T))
    cp.arr[,,1:J] <- aperm(cp.temp, c(1,3,2))
    cp.arr[,,J+1] <- 1 - apply(cp.arr[,,1:J], 1:2, sum, na.rm=TRUE)

    simList <- list()
    for(s in 1:nsim) {
        switch(mixture, 
            P = M <- rpois(n=n, lambda=lam),
            NB = M <- rnbinom(n=n, mu=lam, size=exp(coef(object, type="alpha"))))

        N <- rbinom(n*T, size=M, prob=phi.mat)
        N <- matrix(N, nrow=n, ncol=T, byrow=TRUE)
    
        y <- array(NA, c(n, J, T))
        for(i in 1:n)
            for(t in 1:T)
                y[i,,t] <- drop(rmultinom(1, N[i,t], cp.arr[i,t,]))[1:J] 
        simList[[s]] <- matrix(y, nrow=n, ncol=J*T) # note, byrow=F
        }
    return(simList)
})

          



############################## PARAMETRIC BOOTSTRAP ###########################

setGeneric("parboot",
           def = function(object, ...) {
             standardGeneric("parboot")
           })


setClass("parboot",
         representation(call = "call",
                        t0 = "numeric",
                        t.star = "matrix"))
         

setMethod("parboot", "unmarkedFit", 
    function(object, statistic=SSE, nsim=10, report, ...) 
{
    statistic <- match.fun(statistic)
    call <- match.call(call = sys.call(-1))
    formula <- object@formula
    umf <- getData(object)
    y <- getY(umf)
    if(class(object) %in% c("unmarkedFitOccu", "unmarkedFitOccuRN", 
        "unmarkedFitColExt"))
            y <- truncateToBinary(y)
    ests <- as.numeric(coef(object))
    t0 <- statistic(object, ...)
    lt0 <- length(t0)
    t.star <- matrix(NA, nsim, lt0)
    if(!is.null(names(t0)))
        colnames(t.star) <- names(t0)
    else colnames(t.star) <- paste("t*", 1:lt0, sep="")
    if(!missing(report))
        cat("t0 =", t0, "\n")      
    fits <- list()
    simdata <- umf
    simList <- simulate(object, nsim = nsim, na.rm = FALSE)
    for(i in 1:nsim) {
        y.sim <- simList[[i]]
        is.na(y.sim) <- is.na(y)
        simdata@y <- y.sim
        fits[[i]] <- update(object, data=simdata, starts=ests, se=FALSE, ...)
        t.star[i,] <- statistic(fits[[i]], ...)
        if(!missing(report)) {
            if(nsim > report && i %in% seq(report, nsim, by=report))
                cat(paste(round(t.star[(i-(report-1)):i,], 1), collapse=", "), 
                    fill=TRUE)
            flush.console()
            }
        }
    out <- new("parboot", call=call, t0 = t0, t.star = t.star)
    return(out)
})






setMethod("show", "parboot", function(object) 
{
    t.star <- object@t.star
    t0 <- object@t0
    nsim <- nrow(t.star)
    biasMat <- pMat <- matrix(NA, nsim, length(t0))
    for(i in 1:nsim) {
        biasMat[i,] <- t0 - t.star[i,]
        pMat[i,] <- abs(t.star[i,] - 1) > abs(t0 - 1)
        }
    bias <- colMeans(biasMat)
    bias.se <- apply(biasMat, 2, sd)
    p.val <- colSums(pMat) / (1 + nsim)
    stats <- data.frame("t0" = t0, "mean(t0 - t_B)" = bias, 
        "StdDev(t0 - t_B)" = bias.se, "Pr(t_B > t0)" = p.val, 
        check.names = FALSE)
    cat("\nCall:", deparse(object@call, width=500), fill=T)
    cat("\nParametric Bootstrap Statistics:\n")
    print(stats, digits=3)
    cat("\nt_B quantiles:\n")
    print(t(apply(t.star, 2, quantile, 
        probs=c(0, 2.5, 25, 50, 75, 97.5, 100) / 100)), digits=2)
    cat("\nt0 = Original statistic compuated from data\n")
    cat("t_B = Vector of bootstrap samples\n\n")
})




setMethod("plot", signature(x="parboot", y="missing"), 
    function(x, y, xlab, main = "Parametric Bootstrapped Samples", xlim, ...)
{
    t.star <- x@t.star
    t0 <- x@t0
    for(i in 1:length(t0)) {
        if(missing(xlab))
            xlab <- colnames(t.star)[i]
        h <- hist(t.star[,i], plot = FALSE)
        if(missing(xlim))
            xl <- c(min(h$breaks[1], t0[i]), max(max(h$breaks), t0[i]))
        hist(t.star[,i], xlab=xlab, xlim = xl, main = main, ...)
        abline(v=t0[i], lty=2)
        devAskNewPage(ask = TRUE)
        }
    devAskNewPage(ask = FALSE)
})


# ----------------------- Nonparametric bootstrapping --------------------------

## nonparboot return entire list of fits... 
##  they will be processed by vcov, confint, etc.

setGeneric("nonparboot", 
    function(object, B = 0, ...) {standardGeneric("nonparboot")})


setMethod("nonparboot", "unmarkedFit",
          function(object, B = 0, keepOldSamples = TRUE, bsType, ...) {
            bsType <- match.arg(bsType, c("site", "both"))
            if (identical(B, 0) && !is.null(object@bootstrapSamples)) {
              return(object)
            }
            if (B <= 0 && is.null(object@bootstrapSamples)) {
              stop("B must be greater than 0 when fit has no bootstrap samples.")
            }
            data <- object@data
            formula <- object@formula
            designMats <- getDesign(data, formula)  # bootstrap only after removing sites
            removed.sites <- designMats$removed.sites
            data <- data[-removed.sites,]
            y <- getY(data)
            colnames(y) <- NULL
            data@y <- y
            M <- numSites(data)
            boot.iter <- function() {
              sites <- sort(sample(1:M, M, replace = TRUE))
              data.b <- data[sites,]
              y <- getY(data.b)
              if (bsType == "both") {
                obs.per.site <- alply(y, 1, function(row) {
                  which(!is.na(row))
                })
                obs <- lapply(obs.per.site, 
                    function(obs) sample(obs, replace = TRUE))
                data.b <- data.b[obs]
              }
              fm <- update(object, data = data.b, se = FALSE)
              return(fm)
            }
            if (!keepOldSamples) {
              object@bootstrapSamples <- NULL
            }               
            object@bootstrapSamples <- c(object@bootstrapSamples,
                replicate(B, boot.iter(), simplify = FALSE))
            coefs <- t(sapply(object@bootstrapSamples, function(x) coef(x)))
            v <- cov(coefs)
            object@covMatBS <- v
            inds <- .estimateInds(object)
            for (est in names(inds)) {
              v.est <- v[inds[[est]], inds[[est]], drop = FALSE]
              object@estimates@estimates[[est]]@covMatBS <- v.est
            }
            object
          })

setMethod("nonparboot", "unmarkedFitOccu",
    function(object, B = 0, keepOldSamples = TRUE, ...) 
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples, bsType="both")
})


setMethod("nonparboot", "unmarkedFitPCount",
    function(object, B = 0, keepOldSamples = TRUE, ...) 
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples, bsType="both")
})


setMethod("nonparboot", "unmarkedFitMPois",
    function(object, B = 0, keepOldSamples = TRUE, ...) 
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples, bsType="site")
})


setMethod("nonparboot", "unmarkedFitDS",
    function(object, B = 0, keepOldSamples = TRUE, ...) 
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples, bsType="site")
})


setMethod("nonparboot", "unmarkedFitOccuRN",
    function(object, B = 0, keepOldSamples = TRUE, ...) 
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples, bsType="both")
})


setMethod("nonparboot", "unmarkedFitColExt",
          function(object, B = 0, keepOldSamples = TRUE, ...) {
            if (identical(B, 0) && !is.null(object@bootstrapSamples)) {
              return(object)
            }
            if (B <= 0 && is.null(object@bootstrapSamples)) {
              stop("B must be greater than 0 when fit has no bootstrap samples.")
            }
            data <- object@data
            psiParms <- coef(object, 'psi')
            detParms <- coef(object, 'det')
            colParms <- coef(object, 'col')
            extParms <- coef(object, 'ext')
          
            # bootstrap only after removing sites
            designMats <- getDesign(object@data, formula = object@formula)   
            removed.sites <- designMats$removed.sites
            data <- data[-removed.sites,]
            y <- getY(data)
            colnames(y) <- NULL
            data@y <- y
            M <- numSites(data)
            boot.iter <- function() {
              sites <- sort(sample(1:M, M, replace = TRUE))
              data.b <- data[sites,]
              y <- getY(data.b)
              fm <- update(object, data = data.b, se = FALSE)
              return(fm)
            }
            if (!keepOldSamples) {
              object@bootstrapSamples <- NULL
            }               
            object@bootstrapSamples <- c(object@bootstrapSamples,
                replicate(B, boot.iter(), simplify = FALSE))
            coefs <- t(sapply(object@bootstrapSamples, function(x) coef(x)))
            v <- cov(coefs)
            object@covMatBS <- v
            inds <- .estimateInds(object)
            for (est in names(inds)) {
              v.est <- v[inds[[est]], inds[[est]], drop = FALSE]
              object@estimates@estimates[[est]]@covMatBS <- v.est
            }
            smoothed.occ <- t(sapply(object@bootstrapSamples, 
                function(x) x@smoothed.mean[1,]))
            smoothed.unocc <- t(sapply(object@bootstrapSamples, 
                function(x) x@smoothed.mean[2,]))
            object@smoothed.mean.bsse <- rbind(sqrt(diag(cov(smoothed.occ))),
                                               sqrt(diag(cov(smoothed.unocc))))
            projected.occ <- t(sapply(object@bootstrapSamples, 
                function(x) x@projected.mean[1,]))
            projected.unocc <- t(sapply(object@bootstrapSamples, 
                function(x) x@projected.mean[2,]))
            object@projected.mean.bsse <- rbind(sqrt(diag(cov(projected.occ))),
                                               sqrt(diag(cov(projected.unocc))))
            object
          })

################################# Helper functions #############################

## A helper function to return a list of indices for each estimate type
## 


.estimateInds <- function(umf) {
  ## get length of each estimate
  estimateLengths <- sapply(umf@estimates@estimates, function(est) {
    length(coef(est))
  })
  ## recurse function to generate list of indices
  estimateInds <- function(type) {
    if(type==1) {
      return(list(seq(length=estimateLengths[1])))
    } else {
      prev.list <- estimateInds(type-1)
      prev.max <- max(prev.list[[type-1]])
      return(c(prev.list, list(seq(prev.max+1, prev.max+estimateLengths[type]))))
    }
  }
  retlist <- estimateInds(length(estimateLengths))
  names(retlist) <- names(umf)
  retlist
}
