
#  This function estimates the standard occupancy model of MacKenzie et al (2002).

occu <- function(formula, data, knownOcc = numeric(0), starts, method = "BFGS", 
	control = list(), se = TRUE)
{
	if(!is(data, "unmarkedFrameOccu")) 
    stop("Data is not an unmarkedFrameOccu object.")
		
	designMats <- getDesign(data, formula)
	X <- designMats$X; V <- designMats$V; y <- designMats$y 
  removed <- designMats$removed.sites
  X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
  if (is.null(X.offset)) {
      X.offset <- rep(0, nrow(X))
      }
  if (is.null(V.offset)) {
      V.offset <- rep(0, nrow(V))
      }

	y <- truncateToBinary(y)
	J <- ncol(y)
	M <- nrow(y)

	## convert knownOcc to logical so we can subset correctly to handle NAs.
	knownOccLog <- rep(FALSE, numSites(data))
	knownOccLog[knownOcc] <- TRUE
	knownOccLog <- knownOccLog[-removed]
	
	occParms <- colnames(X)
	detParms <- colnames(V)
	nDP <- ncol(V)
	nOP <- ncol(X)

	nP <- nDP + nOP
	if(!missing(starts) && length(starts) != nP)
	   stop(paste("The number of starting values should be", nP))
	
	yvec <- as.numeric(t(y))
	navec <- is.na(yvec)
	nd <- ifelse(rowSums(y,na.rm=TRUE) == 0, 1, 0) # no det at site i indicator

	nll <- function(params) {
		psi <- plogis(X %*% params[1 : nOP] + X.offset)
		psi[knownOccLog] <- 1
		pvec <- plogis(V %*% params[(nOP + 1) : nP] + V.offset)
		cp <- (pvec^yvec) * ((1 - pvec)^(1 - yvec))
		cp[navec] <- 1  # so that NA's don't modify likelihood
		cpmat <- matrix(cp, M, J, byrow = TRUE) # put back into matrix to multiply appropriately
		loglik <- log(rowProds(cpmat) * psi + nd * (1 - psi))
		-sum(loglik)
	}

	if(missing(starts)) starts <- rep(0, nP)	#rnorm(nP)
	fm <- optim(starts, nll, method = method, control = control, hessian = se)
	opt <- fm
	if(se) {
		tryCatch(covMat <- solve(fm$hessian),
			error=function(x) stop(simpleError("Hessian is singular.  Try using fewer covariates.")))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
	ests <- fm$par
	fmAIC <- 2 * fm$value + 2 * nP #+ 2*nP*(nP + 1)/(M - nP - 1)
	names(ests) <- c(occParms, detParms)

	state <- unmarkedEstimate(name = "Occupancy", short.name = "psi",
		estimates = ests[1:nOP],
		covMat = as.matrix(covMat[1:nOP,1:nOP]), invlink = "logistic",
		invlinkGrad = "logistic.grad")

	det <- unmarkedEstimate(name = "Detection", short.name = "p",
		estimates = ests[(nOP + 1) : nP],
		covMat = as.matrix(covMat[(nOP + 1) : nP, (nOP + 1) : nP]), 
		invlink = "logistic", invlinkGrad = "logistic.grad")

	estimateList <- unmarkedEstimateList(list(state=state, det=det))

	umfit <- new("unmarkedFitOccu", fitType = "occu", call = match.call(), 
		formula = formula, data = data, sitesRemoved = designMats$removed.sites, 
		estimates = estimateList, AIC = fmAIC, opt = opt, negLogLike = fm$value,
		nllFun = nll, knownOcc = knownOccLog)

	return(umfit)
}
