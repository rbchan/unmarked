

setGeneric("getDesign", function(umf, ...) standardGeneric("getDesign"))
setGeneric("handleNA", function(umf, ...) standardGeneric("handleNA"))



# unmarkedFrame

setMethod("getDesign", "unmarkedFrame",
    function(umf, formula, na.rm=TRUE) 
    {
	detformula <- as.formula(formula[[2]])
	stateformula <- as.formula(paste("~", formula[3], sep=""))
	detVars <- all.vars(detformula)
	
	M <- numSites(umf)
	R <- obsNum(umf)
	
	## Compute state design matrix
	if(is.null(siteCovs(umf))) {
		siteCovs <- data.frame(placeHolder = rep(1, M))
	} else {
		siteCovs <- siteCovs(umf)
	}
	X.mf <- model.frame(stateformula, siteCovs, na.action = NULL)
	X <- model.matrix(stateformula, X.mf)
        X.offset <- as.vector(model.offset(X.mf))
        if (!is.null(X.offset)) {
          X.offset[is.na(X.offset)] <- 0
        }

	## Compute detection design matrix
	if(is.null(obsCovs(umf))) {
		obsCovs <- data.frame(placeHolder = rep(1, M*R))
	} else {
		obsCovs <- obsCovs(umf)
	}
	
	## Record future column names for obsCovs
	colNames <- c(colnames(obsCovs), colnames(siteCovs))
	
	## add site Covariates at observation-level
	obsCovs <- cbind(obsCovs, siteCovs[rep(1:M, each = R),])
	colnames(obsCovs) <- colNames
	
	## add observation number if not present
	if(!("obsNum" %in% names(obsCovs))) {
		obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
	}
	
	V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
	V <- model.matrix(detformula, V.mf)
        V.offset <- as.vector(model.offset(V.mf))
	if (!is.null(V.offset)) {
          V.offset[is.na(V.offset)] <- 0
        }

	if (na.rm) {
          out <- handleNA(umf, X, X.offset, V, V.offset)
          y <- out$y
          X <- out$X
          X.offset <- out$X.offset
          V <- out$V
          V.offset <- out$V.offset
          removed.sites <- out$removed.sites
        } else {
          y=getY(umf)
          removed.sites=integer(0)
        }

	return(list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset, 
		removed.sites = removed.sites))
	})
	
	
	

setMethod("handleNA", "unmarkedFrame", function(umf, X, X.offset, V, V.offset) 
    {
	obsToY <- obsToY(umf)
	if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")
	
	J <- numY(umf)
	R <- obsNum(umf)
	M <- numSites(umf)
	
	X.long <- X[rep(1:M, each = J),]
	X.long.na <- is.na(X.long)
	
	V.long.na <- apply(V, 2, function(x) {
				x.mat <- matrix(x, M, R, byrow = TRUE)
				x.mat <- is.na(x.mat)
				x.mat <- x.mat %*% obsToY
				x.long <- as.vector(t(x.mat))
				x.long == 1
			})
	V.long.na <- apply(V.long.na, 1, any)
	
	y.long <- as.vector(t(getY(umf)))
	y.long.na <- is.na(y.long)
	
	covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)
	
	## are any NA in covs not in y already?
	y.new.na <- covs.na & !y.long.na
	
	if(sum(y.new.na) > 0) {
		y.long[y.new.na] <- NA
		warning("Some observations have been discarded because corresponding covariates were missing.", call. = FALSE)
	}
	
	y <- matrix(y.long, M, J, byrow = TRUE)
	sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))
	
	num.to.remove <- sum(sites.to.remove)
	if(num.to.remove > 0) {
		y <- y[!sites.to.remove, ,drop = FALSE]
		X <- X[!sites.to.remove, ,drop = FALSE]
                X.offset <- X.offset[!sites.to.remove]
		V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
                V.offset <- V.offset[!sites.to.remove[rep(1:M, each = R)], ]
		warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
	}
	
	list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
             removed.sites = which(sites.to.remove))
    })
    
    


# UnmarkedMultFrame    
	



setMethod("getDesign", "unmarkedMultFrame", 
    function(umf, formula, na.rm = TRUE) {

    aschar1 <- as.character(formula)
    aschar2 <- as.character(formula[[2]])
    aschar3 <- as.character(formula[[2]][[2]])    

    detformula <- as.formula(paste(aschar1[1], aschar1[3]))
    epsformula <- as.formula(paste(aschar2[1], aschar2[3]))
    gamformula <- as.formula(paste(aschar3[1], aschar3[3]))
    psiformula <- as.formula(formula[[2]][[2]][[2]])
    
    detVars <- all.vars(detformula)
  
    M <- numSites(umf)
    R <- obsNum(umf)
    nY <- umf@numPrimary
  
    ## Compute phi design matrices
    if(is.null(umf@yearlySiteCovs)) {
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*nY))
    } else {
        yearlySiteCovs <- umf@yearlySiteCovs
    }
    ## in order to drop factor levels that only appear in last year,
    ## replace last year with NAs and use drop=TRUE
    yearlySiteCovs[seq(nY,M*nY,by=nY),] <- NA
    yearlySiteCovs <- as.data.frame(lapply(yearlySiteCovs, function(x) {
        x[,drop = TRUE]
        }))
    ## add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = nY),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }
    X.mf.gam <- model.frame(gamformula, yearlySiteCovs, na.action = NULL)
    X.gam <- model.matrix(gamformula, X.mf.gam)
    X.mf.eps <- model.frame(epsformula, yearlySiteCovs, na.action = NULL)
    X.eps <- model.matrix(epsformula, X.mf.eps)
  
    ## Compute site-level design matrix for psi
    if(is.null(siteCovs(umf))) {
        siteCovs <- data.frame(placeHolder = rep(1, M))
    } else {
        siteCovs <- siteCovs(umf)
    }
    W.mf <- model.frame(psiformula, siteCovs, na.action = NULL)
    W <- model.matrix(psiformula, W.mf)

    #  ## impute missing yearlySiteCovs across years as average
    #  X <- t(apply(X, 1, function(x) {
    #            out <- x
    #            out[is.na(x)] <- mean(x)
    #          }))
  
	## Compute detection design matrix
	if(is.null(obsCovs(umf))) {
		obsCovs <- data.frame(placeHolder = rep(1, M*R))
	} else {
		obsCovs <- obsCovs(umf)
	}
	
	## add site and yearlysite covariates at observation-level
	## DOUBLE CHECK THIS. Shouldn't it be:
  ## yearlySiteCovs[rep(1:(M*nY), each = J),] where J is R/nY
	obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*nY), each = R),],
                         siteCovs[rep(1:M, each = R), ])
	
	## add observation number if not present
	if(!("obsNum" %in% names(obsCovs))) {
		obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
	}
	
	V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
	V <- model.matrix(detformula, V.mf)
	
	if(na.rm)
		out <- handleNA(umf, X.gam, X.eps, W, V)
	else
		out <- list(y=getY(umf), X.gam=X.gam, X.eps=X.eps,
                            W=W,V=V,
				removed.sites=integer(0))
	
	return(list(y = out$y, X.eps = out$X.eps, X.gam = out$X.gam,
                    W = out$W, V = out$V,
                    removed.sites = out$removed.sites))
})






setMethod("handleNA", "unmarkedMultFrame", function(umf, X.gam, X.eps, W, V) 
    {
	obsToY <- obsToY(umf)
	if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")
	
	R <- obsNum(umf)
	M <- numSites(umf)
	nY <- umf@numPrimary
	J <- numY(umf) / nY
	
    ## treat both X's and W together
    X <- cbind(X.gam, X.eps, W[rep(1:M, each = nY), ])

	X.na <- is.na(X)
	X.na[seq(nY,M*nY,by=nY),] <- FALSE  ## final years are unimportant (not used).
	X.long.na <- X.na[rep(1:(M*nY), each = J),]
	
	V.long.na <- apply(V, 2, function(x) {
				x.mat <- matrix(x, M, R, byrow = TRUE)
				x.mat <- is.na(x.mat)
				x.mat <- x.mat %*% obsToY
				x.long <- as.vector(t(x.mat))
				x.long == 1
			})
	V.long.na <- apply(V.long.na, 1, any)
	
	y.long <- as.vector(t(getY(umf)))
	y.long.na <- is.na(y.long)
	
	covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)
	
	## are any NA in covs not in y already?
	y.new.na <- covs.na & !y.long.na
	
	if(sum(y.new.na) > 0) {
		y.long[y.new.na] <- NA
		warning("Some observations have been discarded because correspoding covariates were missing.", call. = FALSE)
	}
	
	y <- matrix(y.long, M, numY(umf), byrow = TRUE)
	sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))
	
	num.to.remove <- sum(sites.to.remove)
	if(num.to.remove > 0) {
		y <- y[!sites.to.remove, ,drop = FALSE]
		X.gam <- X.gam[!sites.to.remove[rep(1:M, each = J)], ,drop = FALSE]
                X.eps <- X.eps[!sites.to.remove[rep(1:M, each = J)], ,drop = FALSE]
                W <- X[!sites.to.remove, drop = FALSE]
		V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
		warning(paste(num.to.remove,"sites have been discarded because of missing data."))
	}
	list(y = y, X.gam = X.gam, X.eps = X.eps, W = W,
             V = V,
             removed.sites = which(sites.to.remove))
    })
    
    
    
    
    
    
    
    
    
    
# UnmarkedFrameGMN    
	



setMethod("getDesign", "unmarkedFrameGMM", 
    function(umf, formula, na.rm = TRUE) 
{
    ac1 <- as.character(formula)
    ac2 <- as.character(formula[[2]])

    detformula <- as.formula(paste(ac1[1], ac1[3]))
    phiformula <- as.formula(paste(ac2[1], ac2[3]))
    lamformula <- as.formula(formula[[2]][[2]])
    
    detVars <- all.vars(detformula)
  
    M <- numSites(umf)
    T <- umf@numPrimary
    R <- obsNum(umf)
    J <- R/T
  
    ## Compute phi design matrices
    if(is.null(umf@yearlySiteCovs)) {
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*T))
    } else yearlySiteCovs <- umf@yearlySiteCovs

    # add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = T),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    Xphi.mf <- model.frame(phiformula, yearlySiteCovs, na.action = NULL)
    Xphi <- model.matrix(phiformula, Xphi.mf)
    Xphi.offset <- as.vector(model.offset(Xphi.mf))
    if(!is.null(Xphi.offset)) Xphi.offset[is.na(Xphi.offset)] <- 0
  
    # Compute site-level design matrix for lambda
    if(is.null(siteCovs(umf))) {
        siteCovs <- data.frame(placeHolder = rep(1, M))
    } else siteCovs <- siteCovs(umf)
    Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
    Xlam <- model.matrix(lamformula, Xlam.mf)
    Xlam.offset <- as.vector(model.offset(Xlam.mf))
    if(!is.null(Xlam.offset)) Xlam.offset[is.na(Xlam.offset)] <- 0

    # Compute detection design matrix
    if(is.null(obsCovs(umf))) {
        obsCovs <- data.frame(placeHolder = rep(1, M*R))
    } else obsCovs <- obsCovs(umf)
	
    # add site and yearlysite covariates at observation-level
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*T), each = J),],
        siteCovs[rep(1:M, each = R), ])
	
    # add observation number if not present
    # FIXME. obs is often used as abbreviation for observer. 
    # This can cause conflict
    if(!("obsNum" %in% names(obsCovs)))
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))
	
    Xdet.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    Xdet <- model.matrix(detformula, Xdet.mf)
    Xdet.offset <- as.vector(model.offset(Xdet.mf))
    if(!is.null(Xdet.offset)) Xdet.offset[is.na(Xdet.offset)] <- 0
	
    if(na.rm)
        out <- handleNA(umf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet, 
            Xdet.offset)
    else
        out <- list(y=getY(umf), Xlam=Xlam, Xlam.offset = Xlam.offset, 
            Xphi=Xphi, Xphi.offset = Xphi.offset, Xdet=Xdet,
				    removed.sites=integer(0))
	
    return(list(y = out$y, Xlam = out$Xlam, Xphi = out$Xphi, Xdet = out$Xdet, 
        Xlam.offset = out$Xlam.offset, Xphi.offset = out$Xphi.offset, 
        Xdet.offset = out$Xdet.offset, removed.sites = out$removed.sites))
})






setMethod("handleNA", "unmarkedFrameGMM", 
    function(umf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet, Xdet.offset) 
{

    obsToY <- obsToY(umf)
        if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

    M <- numSites(umf)
    T <- umf@numPrimary
    R <- obsNum(umf)
    J <- R/T            # numY(umf) / T
	
    # treat Xphi and Xlam together
    X <- cbind(Xphi, Xlam[rep(1:M, each = T), ])

    X.na <- is.na(X)
    X.long.na <- X.na[rep(1:(M*T), each = J),]
	
    Xdet.long.na <- apply(Xdet, 2, function(x) {
        x.mat <- matrix(x, M, R, byrow = TRUE)
        x.mat <- is.na(x.mat)
        x.mat <- x.mat %*% obsToY
        x.long <- as.vector(t(x.mat))
        x.long == 1
        })

    Xdet.long.na <- apply(Xdet.long.na, 1, any)
	
    y.long <- as.vector(t(getY(umf)))
    y.long.na <- is.na(y.long)
	
    covs.na <- apply(cbind(X.long.na, Xdet.long.na), 1, any)
	
    ## are any NA in covs not in y already?
    y.new.na <- covs.na & !y.long.na
	
    if(sum(y.new.na) > 0) {
    		y.long[y.new.na] <- NA
    		warning("Some observations have been discarded because correspoding covariates were missing.", call. = FALSE)
	     }
	
    y <- matrix(y.long, M, numY(umf), byrow = TRUE)
    sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))
	
    num.to.remove <- sum(sites.to.remove)
    if(num.to.remove > 0) {
    		y <- y[!sites.to.remove, ,drop = FALSE]
        Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
    		Xlam.offset <- Xlam.offset[!sites.to.remove]
    		Xphi <- Xphi[!sites.to.remove[rep(1:M, each = J)], ,drop = FALSE]
    		Xphi.offset <- Xphi.offset[!sites.to.remove[rep(1:M, each = J)]]
        Xdet <- Xdet[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
        Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=R)]]
    		warning(paste(num.to.remove, 
            "sites have been discarded because of missing data."), call.=FALSE)
	     }
    list(y = y, Xlam = Xlam, Xlam.offset = Xlam.offset, Xphi = Xphi, 
        Xphi.offset = Xphi.offset, Xdet = Xdet, Xdet.offset = Xdet.offset,
        removed.sites = which(sites.to.remove))
})
    
