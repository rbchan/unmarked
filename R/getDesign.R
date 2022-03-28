

setGeneric("getDesign", function(umf, ...) standardGeneric("getDesign"))
setGeneric("handleNA", function(umf, ...) standardGeneric("handleNA"))



# unmarkedFrame

setMethod("getDesign", "unmarkedFrame",
    function(umf, formula, na.rm=TRUE)
{
    detformula <- lme4::nobars(as.formula(formula[[2]]))
    stateformula <- lme4::nobars(as.formula(paste("~", formula[3], sep="")))
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

    #Generate random effects indicator matrices
    sf_rand <- as.formula(paste("~", formula[3], sep=""))
    df_rand <- as.formula(formula[[2]])
    Z_state <- get_Z(sf_rand, siteCovs)
    Z_det <- get_Z(df_rand, obsCovs)

    if (na.rm) {
        out <- handleNA(umf, X, X.offset, V, V.offset, Z_state, Z_det)
        y <- out$y
        X <- out$X
        X.offset <- out$X.offset
        V <- out$V
        V.offset <- out$V.offset
        Z_state <- out$Z_state
        Z_det <- out$Z_det
        removed.sites <- out$removed.sites
    } else {
        y=getY(umf)
        removed.sites=integer(0)
    }

    if (is.null(X.offset)) X.offset <- rep(0, nrow(X))
    if (is.null(V.offset)) V.offset <- rep(0, nrow(V))

    return(list(y = y, X = X, X.offset = X.offset, V = V,
                V.offset = V.offset,
                Z_state = Z_state, Z_det = Z_det,
                removed.sites = removed.sites))
})


setMethod("handleNA", "unmarkedFrame",
          function(umf, X, X.offset, V, V.offset, Z_state, Z_det)
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
        x.long > 0
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
        if(nrow(Z_state)>0) Z_state <- Z_state[!sites.to.remove, , drop=FALSE]
        if(nrow(Z_det)>0) Z_det <- Z_det[!sites.to.remove[rep(1:M, each = R)],,drop=FALSE]
        warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
        }

    list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
         Z_state = Z_state, Z_det = Z_det,
        removed.sites = which(sites.to.remove))
    })


# unmarkedFrameOccuFP
# like the occuFP but there are 3 osbservation formula which are stored in V (true positive detections),
# U (false positive detections), and W (b or probability detetion is certain)


setMethod("getDesign", "unmarkedFrameOccuFP",
          function(umf, detformula,FPformula,Bformula = ~.,stateformula, na.rm=TRUE)
          {

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


            U.mf <- model.frame(FPformula, obsCovs, na.action = NULL)
            U <- model.matrix(FPformula, U.mf)
            U.offset <- as.vector(model.offset(U.mf))
            if (!is.null(U.offset)) {
              U.offset[is.na(U.offset)] <- 0
            }

            W.mf <- model.frame(Bformula, obsCovs, na.action = NULL)
            W <- model.matrix(Bformula, W.mf)
            W.offset <- as.vector(model.offset(W.mf))
            if (!is.null(W.offset)) {
              W.offset[is.na(W.offset)] <- 0
            }

            if (na.rm) {
              out <- handleNA(umf, X, X.offset, V,V.offset, U, U.offset, W, W.offset)
              y <- out$y
              X <- out$X
              X.offset <- out$X.offset
              V <- out$V
              V.offset <- out$V.offset
              U <- out$U
              U.offset <- out$U.offset
              U <- out$U
              U.offset <- out$U.offset
              removed.sites <- out$removed.sites
            } else {
              y=getY(umf)
              removed.sites=integer(0)
            }



            return(list(y = y, X = X, X.offset = X.offset, V = V,
                        V.offset = V.offset,U = U, U.offset = U.offset,W = W,
                        W.offset = W.offset, removed.sites = removed.sites))
          })


setMethod("handleNA", "unmarkedFrameOccuFP",
          function(umf, X, X.offset, V, V.offset, U, U.offset, W, W.offset)
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
              x.long > 0
            })
            V.long.na <- apply(V.long.na, 1, any)

            U.long.na <- apply(U, 2, function(x) {
              x.mat <- matrix(x, M, R, byrow = TRUE)
              x.mat <- is.na(x.mat)
              x.mat <- x.mat %*% obsToY
              x.long <- as.vector(t(x.mat))
              x.long > 0
            })
            U.long.na <- apply(U.long.na, 1, any)

            W.long.na <- apply(W, 2, function(x) {
              x.mat <- matrix(x, M, R, byrow = TRUE)
              x.mat <- is.na(x.mat)
              x.mat <- x.mat %*% obsToY
              x.long <- as.vector(t(x.mat))
              x.long > 0
            })
            W.long.na <- apply(W.long.na, 1, any)

            y.long <- as.vector(t(getY(umf)))
            y.long.na <- is.na(y.long)

            covs.na <- apply(cbind(X.long.na, V.long.na, U.long.na, W.long.na), 1, any)

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
              U <- U[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
              U.offset <- U.offset[!sites.to.remove[rep(1:M, each = R)], ]
              W <- W[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
              W.offset <- W.offset[!sites.to.remove[rep(1:M, each = R)], ]
              warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
            }

            list(y = y, X = X, X.offset = X.offset, V = V, V.offset = V.offset,
                 U = U, U.offset = U.offset, W = W, W.offset = W.offset,
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
    J <- R / nY

    ## Compute phi design matrices
    if(is.null(umf@yearlySiteCovs)) {
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*nY))
    } else {
        yearlySiteCovs <- umf@yearlySiteCovs
    }
    ## add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = nY),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    ## Compute site-level design matrix for psi
    if(is.null(siteCovs(umf)))
        siteCovs <- data.frame(placeHolder = rep(1, M))
    else
        siteCovs <- siteCovs(umf)

    W.mf <- model.frame(psiformula, siteCovs, na.action = NULL)
    if(!is.null(model.offset(W.mf)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    W <- model.matrix(psiformula, W.mf)


    ## Compute detection design matrix
    if(is.null(obsCovs(umf)))
        obsCovs <- data.frame(placeHolder = rep(1, M*R))
    else
        obsCovs <- obsCovs(umf)

    ## add site and yearlysite covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(yearlySiteCovs))
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*nY), each = J),])
    colnames(obsCovs) <- cnames

    ## add observation number if not present
    if(!("obsNum" %in% names(obsCovs)))
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:R, M)))

    V.mf <- model.frame(detformula, obsCovs, na.action = NULL)
    if(!is.null(model.offset(V.mf)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    V <- model.matrix(detformula, V.mf)

    ## in order to drop factor levels that only appear in last year,
    ## replace last year with NAs and use drop=TRUE
    yearlySiteCovs[seq(nY,M*nY,by=nY),] <- NA
    yearlySiteCovs <- as.data.frame(lapply(yearlySiteCovs, function(x) {
        x[,drop = TRUE]
        }))

    X.mf.gam <- model.frame(gamformula, yearlySiteCovs, na.action = NULL)
    if(!is.null(model.offset(X.mf.gam)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    X.gam <- model.matrix(gamformula, X.mf.gam)
    X.mf.eps <- model.frame(epsformula, yearlySiteCovs, na.action = NULL)
    if(!is.null(model.offset(X.mf.eps)))
        stop("offsets not currently allowed in colext", call.=FALSE)
    X.eps <- model.matrix(epsformula, X.mf.eps)

    if(na.rm)
        out <- handleNA(umf, X.gam, X.eps, W, V)
    else
        out <- list(y=getY(umf), X.gam=X.gam, X.eps=X.eps, W=W, V=V,
            removed.sites=integer(0))

    return(list(y = out$y, X.eps = out$X.eps, X.gam = out$X.gam, W = out$W,
        V = out$V, removed.sites = out$removed.sites))
})






setMethod("handleNA", "unmarkedMultFrame",
    function(umf, X.gam, X.eps, W, V)
{
    obsToY <- obsToY(umf)
    if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

    R <- obsNum(umf)
    M <- numSites(umf)
    nY <- umf@numPrimary
    J <- numY(umf) / nY

    ## treat both X's #######no: and W together
#    X <- cbind(X.gam, X.eps, W[rep(1:M, each = nY), ])
    X <- cbind(X.gam, X.eps)

    X.na <- is.na(X)
    X.na[seq(nY,M*nY,by=nY),] <- FALSE  ## final years are unimportant.
                                        ## not true for W covs!!!
    W.expand <- W[rep(1:M, each=nY),,drop=FALSE]
    W.na <- is.na(W.expand)
    X.na <- cbind(X.na, W.na) # NAs in siteCovs results in removed site

    X.long.na <- X.na[rep(1:(M*nY), each = J),]

    V.long.na <- apply(V, 2, function(x) {
        x.mat <- matrix(x, M, R, byrow = TRUE)
        x.mat <- is.na(x.mat)
        x.mat <- x.mat %*% obsToY
        x.long <- as.vector(t(x.mat))
        x.long > 0
    })
    V.long.na <- apply(V.long.na, 1, any)

    y.long <- as.vector(t(getY(umf)))
    y.long.na <- is.na(y.long)

    # It doesn't make sense to combine X.gam/eps with W here b/c
    # a X.eps does not map correctly to y
    covs.na <- apply(cbind(X.long.na, V.long.na), 1, any)

    ## are any NA in covs not in y already?
    y.new.na <- covs.na & !y.long.na

    if(sum(y.new.na) > 0) {
        y.long[y.new.na] <- NA
        warning("Some observations have been discarded because correspoding covariates were missing.", call. = FALSE)
        }

    y <- matrix(y.long, M, numY(umf), byrow = TRUE)
    sites.to.remove <- apply(y, 1, function(x) all(is.na(x)))
#    Perhaps we need to remove sites that have no data in T=1
#    noDataT1 <- apply(is.na(y[,1:J]), 1, all) #
#    sites.to.remove <- sites.to.remove | noDataT1

    num.to.remove <- sum(sites.to.remove)
    if(num.to.remove > 0) {
        y <- y[!sites.to.remove, ,drop = FALSE]
        X.gam <- X.gam[!sites.to.remove[rep(1:M, each = nY)],,drop = FALSE]
        X.eps <- X.eps[!sites.to.remove[rep(1:M, each = nY)],,drop = FALSE]
        W <- W[!sites.to.remove,, drop = FALSE] # !!! Recent bug fix
        V <- V[!sites.to.remove[rep(1:M, each = R)], ,drop = FALSE]
        warning(paste(num.to.remove,"sites have been discarded because of missing data."), call. = FALSE)
    }
    list(y = y, X.gam = X.gam, X.eps = X.eps, W = W, V = V,
        removed.sites = which(sites.to.remove))
})

# occuMulti

setMethod("getDesign", "unmarkedFrameOccuMulti",
    function(umf, detformulas, stateformulas, maxOrder, na.rm=TRUE, warn=FALSE,
             newdata=NULL, type="state")
{

  #Format formulas
  #Workaround for parameters fixed at 0
  fixed0 <- stateformulas %in% c("~0","0")
  stateformulas[fixed0] <- "~1"

  stateformulas <- lapply(stateformulas,as.formula)
  detformulas <- lapply(detformulas,as.formula)

  #Generate some indices
  S <- length(umf@ylist) # of species
  if(missing(maxOrder)){
    maxOrder <- S
  }
  z <- expand.grid(rep(list(1:0),S))[,S:1] # z matrix
  colnames(z) <- names(umf@ylist)
  M <- nrow(z) # of possible z states

  # f design matrix
  if(maxOrder == 1){
    dmF <- as.matrix(z)
  } else {
    dmF <- model.matrix(as.formula(paste0("~.^",maxOrder,"-1")),z)
  }
  nF <- ncol(dmF) # of f parameters

  J <- ncol(umf@ylist[[1]]) # max # of samples at a site
  N <- nrow(umf@ylist[[1]]) # of sites

  #Check formulas
  if(length(stateformulas) != nF)
    stop(paste(nF,"formulas are required in stateformulas list"))
  if(length(detformulas) != S)
    stop(paste(S,"formulas are required in detformulas list"))

  if(is.null(siteCovs(umf))) {
    site_covs <- data.frame(placeHolderSite = rep(1, N))
  } else {
    site_covs <- siteCovs(umf)
  }

  if(is.null(obsCovs(umf))) {
    obs_covs <- data.frame(placeHolderObs = rep(1, J*N))
  } else {
    obs_covs <- obsCovs(umf)
  }

  #Add site covs to obs covs if we aren't predicting with newdata
  # Record future column names for obsCovs
  col_names <- c(colnames(obs_covs), colnames(site_covs))

  # add site covariates at observation-level
  obs_covs <- cbind(obs_covs, site_covs[rep(1:N, each = J),])
  colnames(obs_covs) <- col_names

  #Re-format ylist
  index <- 1
  ylong <- lapply(umf@ylist, function(x) {
                   colnames(x) <- 1:J
                   x <- cbind(x,site=1:N,species=index)
                   index <<- index+1
                   x
          })
  ylong <- as.data.frame(do.call(rbind,ylong))
  ylong <- reshape(ylong, idvar=c("site", "species"), varying=list(1:J),
                   v.names="value", direction="long")
  ylong <- reshape(ylong, idvar=c("site","time"), v.names="value",
                    timevar="species", direction="wide")
  ylong <- ylong[order(ylong$site, ylong$time), ]

  #Remove missing values
  if(na.rm){
    naSiteCovs <- which(apply(site_covs, 1, function(x) any(is.na(x))))
    if(length(naSiteCovs>0)){
      stop(paste("Missing site covariates at sites:",
                 paste(naSiteCovs,collapse=", ")))
    }

    naY <- apply(ylong, 1, function(x) any(is.na(x)))
    naCov <- apply(obs_covs, 1, function(x) any(is.na(x)))
    navec <- naY | naCov

    sites_with_missingY <- unique(ylong$site[naY])
    sites_with_missingCov <- unique(ylong$site[naCov])

    ylong <- ylong[!navec,,drop=FALSE]
    obs_covs <- obs_covs[!navec,,drop=FALSE]

    no_data_sites <- which(! 1:N %in% ylong$site)
    if(length(no_data_sites>0)){
      stop(paste("All detections and/or detection covariates are missing at sites:",
                  paste(no_data_sites,collapse=", ")))
    }

    if(sum(naY)>0&warn){
      warning(paste("Missing detections at sites:",
                    paste(sites_with_missingY,collapse=", ")))
    }
    if(sum(naCov)>0&warn){
      warning(paste("Missing detection covariate values at sites:",
                    paste(sites_with_missingCov,collapse=", ")))
    }

  }

  #Start-stop indices for sites
  yStart <- c(1,1+which(diff(ylong$site)!=0))
  yStop <- c(yStart[2:length(yStart)]-1,nrow(ylong))

  y <- as.matrix(ylong[,3:ncol(ylong)])

  #Indicator matrix for no detections at a site
  Iy0 <- do.call(cbind, lapply(umf@ylist,
                               function(x) as.numeric(rowSums(x, na.rm=T)==0)))

  #Save formatted covariate frames for use in model frames
  #For predicting with formulas etc
  site_ref <- site_covs
  obs_ref <- obs_covs

  #Assign newdata as the covariate frame if it is provided
  if(!is.null(newdata)){
    if(type == "state"){
      site_covs <- newdata
    } else if(type == "det"){
      obs_covs <- newdata
    }
  }

  #Design matrices + parameter counts
  #For f/occupancy
  fInd <- c()
  sf_no0 <- stateformulas[!fixed0]
  var_names <- colnames(dmF)[!fixed0]
  dmOcc <- lapply(seq_along(sf_no0),function(i){
                    fac_col <- site_ref[, sapply(site_ref, is.factor), drop=FALSE]
                    mf <- model.frame(sf_no0[[i]], site_ref)
                    xlevs <- lapply(fac_col, levels)
                    xlevs <- xlevs[names(xlevs) %in% names(mf)]
                    out <- model.matrix(sf_no0[[i]],
                                        model.frame(stats::terms(mf), site_covs, na.action=stats::na.pass, xlev=xlevs))
                    colnames(out) <- paste('[',var_names[i],'] ',
                                           colnames(out), sep='')
                    fInd <<- c(fInd,rep(i,ncol(out)))
                    out
          })
  fStart <- c(1,1+which(diff(fInd)!=0))
  fStop <- c(fStart[2:length(fStart)]-1,length(fInd))
  occParams <- unlist(lapply(dmOcc,colnames))
  nOP <- length(occParams)

  #For detection
  dInd <- c()
  dmDet <- lapply(seq_along(detformulas),function(i){
                    fac_col <- obs_ref[, sapply(obs_ref, is.factor), drop=FALSE]
                    mf <- model.frame(detformulas[[i]], obs_ref)
                    xlevs <- lapply(fac_col, levels)
                    xlevs <- xlevs[names(xlevs) %in% names(mf)]
                    out <- model.matrix(detformulas[[i]],
                                        model.frame(stats::terms(mf), obs_covs, na.action=stats::na.pass, xlev=xlevs))
                    colnames(out) <- paste('[',names(umf@ylist)[i],'] ',
                                           colnames(out),sep='')
                    dInd <<- c(dInd,rep(i,ncol(out)))
                    out
          })
  dStart <- c(1,1+which(diff(dInd)!=0)) + nOP
  dStop <- c(dStart[2:length(dStart)]-1,length(dInd)+nOP)
  detParams <- unlist(lapply(dmDet,colnames))
  #nD <- length(detParams)

  #Combined
  paramNames <- c(occParams,detParams)
  nP <- length(paramNames)

  mget(c("N","S","J","M","nF","fStart","fStop","fixed0","dmF","dmOcc","dmDet",
         "dStart","dStop","y","yStart","yStop","Iy0","z","nOP","nP","paramNames"))
})

## occuMS

setMethod("getDesign", "unmarkedFrameOccuMS",
    function(umf, psiformulas, phiformulas, detformulas, prm, na.rm=TRUE,
             newdata=NULL, type="psi")
{

  N <- numSites(umf)
  S <- umf@numStates
  T <- umf@numPrimary
  R <- obsNum(umf)
  J <- R / T
  npsi <- S-1 #Number of free psi values
  nphi <- S^2 - S #Number of free phi values
  np <- S * (S-1) / 2 #Number of free p values

  if(length(psiformulas) != npsi){
    stop(paste(npsi,'formulas are required in psiformulas vector'))
  }

  if(is.null(phiformulas)){
    if(prm == 'condbinom') {
      phiformulas <- rep('~1',6)
    } else {
      phiformulas <- rep('~1',S^2-S)
    }
  } else if(T>1){
    if(length(phiformulas)!=nphi){
      stop(paste(nphi,'formulas are required in phiformulas vector. See data@phiOrder for help'))
    }
  }

  if(length(detformulas) != np){
    stop(paste(np,'formulas are required in detformulas vector'))
  }

  #get placeholder for empty covs if necessary
  get_covs <- function(covs, length){
    if(is.null(covs)){
      return(data.frame(placeHolder = rep(1, length)))
    }
    return(covs)
  }

  #Function to create list of design matrices from list of formulas
  get_dm <- function(formulas, covs, namevec, old_covs=NULL){

    ref_covs <- covs
    if(!is.null(old_covs)) ref_covs <- old_covs

    fac_col <- ref_covs[, sapply(ref_covs, is.factor), drop=FALSE]
    xlevs_all <- lapply(fac_col, levels)

    apply_func <- function(i){
      mf <- model.frame(as.formula(formulas[i]), ref_covs)
      xlevs <- xlevs_all[names(xlevs_all) %in% names(mf)]
      out <- model.matrix(as.formula(formulas[i]),
              model.frame(stats::terms(mf), covs,
                          na.action=stats::na.pass, xlev=xlevs))
      #improve these names
      colnames(out) <- paste(namevec[i], colnames(out))
      out
    }

    out <- lapply(seq_along(formulas), apply_func)
    names(out) <- namevec
    out
  }

  #Generate informative names for p
  get_p_names <- function(S, prm){
    if(prm=='condbinom'){
      return(c('p[1]','p[2]','delta'))
    }
    inds <- matrix(NA,nrow=S,ncol=S)
    inds <- lower.tri(inds,diag=TRUE)
    inds[,1] <- FALSE
    inds <- which(inds,arr.ind=TRUE) - 1
    paste0('p[',inds[,2],inds[,1],']')
  }

  #Informative names for phi
  get_phi_names <- function(np, prm){
    if(prm=='condbinom'){
      return(c(paste0('phi[',0:(S-1),']'),paste0('R[',0:(S-1),']')))
    }
    vals <- paste0('phi[',rep(0:(S-1),each=S),rep(0:(S-1),S),']')
    vals <- matrix(vals,nrow=S)
    diag(vals) <- NA
    c(na.omit(as.vector(vals)))
  }

  #Informative names for psi
  get_psi_names <- function(np, prm){
    if(prm=='condbinom'){
      return(c('psi','R'))
    }
    paste0('psi[',1:np,']')
  }

  #Get vector of parameter count indices from a design matrix list
  get_param_inds <- function(dm_list, offset=0){
    apply_func <- function(i){
      rep(i, ncol(dm_list[[i]]))
    }
    ind_vec <- unlist(lapply(seq_along(dm_list), apply_func))
    start_ind <- c(1,1+which(diff(ind_vec)!=0)) + offset
    stop_ind <- c(start_ind[2:length(start_ind)]-1,length(ind_vec)+offset)

    cbind(start_ind, stop_ind)
  }

  #Get param names from dm_list
  get_param_names <- function(dm_list){
    unlist(lapply(dm_list,colnames))
  }

  #Get observations with NAs across design matrices
  get_na_inds <- function(formulas, covs){
    dm_list <- get_dm(formulas, covs, rep("", length(formulas)))
    dm_mat <- Reduce(cbind, dm_list)
    which(apply(dm_mat, 1, function(x) any(is.na(x))))
  }

  site_covs <- get_covs(siteCovs(umf), N)

  y_site_covs <- get_covs(yearlySiteCovs(umf), N*T)

  #Add site covs to yearly site covs
  if(!is.null(umf@siteCovs)){
    y_site_covs <- cbind(y_site_covs,
                         site_covs[rep(1:N, each=T),,drop=FALSE])
  }

  obs_covs <- get_covs(obsCovs(umf), N*R)
  #Add yearly site covs to obs covs
  if(!is.null(umf@siteCovs) | !is.null(umf@yearlySiteCovs)){
    obs_covs <- cbind(obs_covs,
                      y_site_covs[rep(1:(N*T), each=J),,drop=FALSE])
  }

  ## in order to drop factor levels that only appear in last year,
  ## replace last year with NAs and use drop=TRUE
  y_site_covs[seq(T,N*T,by=T),] <- NA
  y_site_covs <- as.data.frame(lapply(y_site_covs, function(x) x[,drop = TRUE]))
  #Actually just remove last year
  y_site_covs <- y_site_covs[-seq(T,N*T,by=T),,drop=FALSE]

  y <- getY(umf)

  #Handle NAs
  removed.sites <- NA
  if(na.rm){

    #Det
    ylong <- as.vector(t(y))
    miss_det_cov <- get_na_inds(detformulas, obs_covs)
    miss_y <- which(is.na(ylong))
    new_na <- miss_det_cov[!miss_det_cov%in%miss_y]
    if(length(new_na)>0){
      warning('Some observations removed because covariates were missing')
      ylong[new_na] <- NA
      y <- matrix(ylong,nrow=N,ncol=R,byrow=T)
    }

    #State
    check_site_na <- function(yrow){
      if(T==1) return(all(is.na(yrow)))
      y_mat <- matrix(yrow, nrow=J)
      pp_na <- apply(y_mat,2,function(x) all(is.na(x)))
      if(any(pp_na)){
        return(TRUE)
      }
      return(FALSE)
    }

    all_y_na <- which(apply(y,1, check_site_na ))
    if(length(all_y_na)>0){
      warning("Some sites removed because all y values in a primary period were missing")
    }
    miss_covs <- get_na_inds(psiformulas, site_covs)
    if(length(miss_covs)>0){
      warning("Some sites removed because site covariates were missing")
    }
    removed.sites <- sort(unique(c(all_y_na,miss_covs)))

    if(T>1){
      ysc_na <- get_na_inds(phiformulas, y_site_covs)
      if(length(ysc_na) > 0){
        stop("Some sites are missing yearly site covs")
      }
    }

    if(length(removed.sites)>0){
      ymap <- as.vector(t(matrix(rep(1:N,each=R),ncol=R,byrow=T)))
      site_covs <- site_covs[-removed.sites,,drop=FALSE]
      obs_covs <- obs_covs[!ymap%in%removed.sites,,drop=FALSE]
      if(T>1){
        ysc_map <- as.vector(t(matrix(rep(1:N,each=(T-1)),ncol=(T-1),byrow=T)))
        y_site_covs <- y_site_covs[!ysc_map%in%removed.sites,,drop=FALSE]
      }
      y <- y[-removed.sites,,drop=FALSE]
      N <- nrow(y)
    }

  }

  site_ref <- site_covs
  ysc_ref <- y_site_covs
  obs_ref <- obs_covs

  #Assign newdata as the covariate frame if it is provided
  if(!is.null(newdata)){
    if(type == "psi"){
      site_covs <- newdata
    } else if(type == "phi"){
      y_site_covs <- newdata
    } else if(type == "det"){
      obs_covs <- newdata
    }
  }

  dm_state <- get_dm(psiformulas, site_covs,
                     get_psi_names(length(psiformulas),prm), site_ref)
  nSP <- length(get_param_names(dm_state))
  state_ind <- get_param_inds(dm_state) #generate ind matrix in function

  nPP <- 0; dm_phi <- list(); phi_ind <- c()
  if(T>1){
    dm_phi <- get_dm(phiformulas, y_site_covs,
                   get_phi_names(length(phiformulas),prm), ysc_ref)
    nPP <- length(get_param_names(dm_phi))
    phi_ind <- get_param_inds(dm_phi, offset=nSP)
  }

  dm_det <- get_dm(detformulas, obs_covs, get_p_names(S,prm), obs_ref)
  det_ind <- get_param_inds(dm_det, offset=(nSP+nPP))

  param_names <- c(get_param_names(dm_state),
                   get_param_names(dm_phi),
                   get_param_names(dm_det))

  mget(c("y","dm_state","state_ind","nSP",
         "dm_phi","phi_ind","nPP",
         "dm_det","det_ind","param_names","removed.sites"))

})

# pcountOpen
#setMethod("getDesign", "unmarkedFramePCOorMMO",
setMethod("getDesign", "unmarkedFramePCO",
    function(umf, formula, na.rm = TRUE)
{
    aschar1 <- as.character(formula)
    aschar2 <- as.character(formula[[2]])
    aschar3 <- as.character(formula[[2]][[2]])
    aschar4 <- as.character(formula[[2]][[2]][[2]])

    iotaformula <- as.formula(paste(aschar1[1], aschar1[3]))
    pformula <- as.formula(paste(aschar2[1], aschar2[3]))
    omformula <- as.formula(paste(aschar3[1], aschar3[3]))
    gamformula <- as.formula(paste(aschar4[1], aschar4[3]))
    lamformula <- as.formula(formula[[2]][[2]][[2]][[2]])

    y <- getY(umf)
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    R <- obsNum(umf) / T
    delta <- umf@primaryPeriod

    if(is.null(umf@yearlySiteCovs))
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*T))
    else
        yearlySiteCovs <- umf@yearlySiteCovs

    ## add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = T),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    if(is.null(siteCovs(umf)))
        siteCovs <- data.frame(placeHolder = rep(1, M))
    else
        siteCovs <- siteCovs(umf)

    Xlam.mf <- model.frame(lamformula, siteCovs, na.action = NULL)
    Xlam <- model.matrix(lamformula, Xlam.mf)
    Xlam.offset <- as.vector(model.offset(Xlam.mf))
    if(!is.null(Xlam.offset))
        Xlam.offset[is.na(Xlam.offset)] <- 0

    if(is.null(obsCovs(umf)))
        obsCovs <- data.frame(placeHolder = rep(1, M*R*T))
    else
        obsCovs <- obsCovs(umf)

    colNames <- c(colnames(obsCovs), colnames(yearlySiteCovs))

    # Add yearlySiteCovs, which contains siteCovs
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*T), each = R),])
    colnames(obsCovs) <- colNames

    if(!("obsNum" %in% names(obsCovs)))
        obsCovs <- cbind(obsCovs, obsNum = as.factor(rep(1:(R*T), M)))

    # Ignore last year of data
    transCovs <- yearlySiteCovs[-seq(T, M*T, by=T),,drop=FALSE]
    for(i in 1:ncol(transCovs))
        if(is.factor(transCovs[,i]))
            transCovs[,i] <- factor(transCovs[,i]) # drop unused levels

    Xiota.mf <- model.frame(iotaformula, transCovs, na.action = NULL)
    Xiota <- model.matrix(iotaformula, Xiota.mf)
    Xiota.offset <- as.vector(model.offset(Xiota.mf))
    if(!is.null(Xiota.offset))
        Xiota.offset[is.na(Xiota.offset)] <- 0
    Xp.mf <- model.frame(pformula, obsCovs, na.action = NULL)
    Xp <- model.matrix(pformula, Xp.mf)
    Xp.offset <- as.vector(model.offset(Xp.mf))
    if(!is.null(Xp.offset))
        Xp.offset[is.na(Xp.offset)] <- 0
    Xgam.mf <- model.frame(gamformula, transCovs, na.action = NULL)
    Xgam <- model.matrix(gamformula, Xgam.mf)
    Xgam.offset <- as.vector(model.offset(Xgam.mf))
    if(!is.null(Xgam.offset))
        Xgam.offset[is.na(Xgam.offset)] <- 0
    Xom.mf <- model.frame(omformula, transCovs, na.action = NULL)
    Xom <- model.matrix(omformula, Xom.mf)
    Xom.offset <- as.vector(model.offset(Xom.mf))
    if(!is.null(Xom.offset))
        Xom.offset[is.na(Xom.offset)] <- 0

    # determine if gamma, omega, and iota are scalar, vector, or matrix valued
    # Runtime is much faster for scalars and vectors
    Xgo <- cbind(Xgam, Xom, Xiota)
    getGOdims <- function(x) {
        xm <- matrix(x, M, T-1, byrow=TRUE)
#        anyNA <- apply(is.na(xm), 1, any)
#        if(all(anyNA))
#            return("matrix")
#        xm <- xm[!anyNA,] # This is not 100% safe
#        nSites <- nrow(xm)
#        if(all(dim(unique(xm, MARGIN=1)) == c(1, T-1)))
#            return("rowvec")
#        else if(all(dim(unique(xm, MARGIN=2)) == c(nSites, 1)))
#            return("colvec")
#        else return("matrix")
        col.table <- apply(xm, 2, table)
        row.table <- apply(xm, 1, table)
        if(is.vector(col.table) & !is.list(col.table)) {
            return("rowvec")
        } else if(is.vector(row.table) & !is.list(row.table)) {
            return("colvec")
        } else
            return("matrix")
        }
    if(isTRUE(all.equal(gamformula,~1)) & isTRUE(all.equal(omformula, ~1)) &
      isTRUE(all.equal(iotaformula, ~1)))
        go.dims <- "scalar"
    else {
        go.dims.vec <- apply(Xgo, 2, getGOdims)
        if(all(go.dims.vec == "rowvec"))
            go.dims <- "rowvec"
        else if(all(go.dims.vec == "colvec"))
            go.dims <- "matrix" ##"colvec"  ## NOTE: Temporary fix to the problem reported with time-only-varying covariates
        else
            go.dims <- "matrix"
    }

    if(na.rm)
        out <- handleNA(umf, Xlam, Xgam, Xom, Xp, Xiota,
            Xlam.offset, Xgam.offset, Xom.offset, Xp.offset, Xiota.offset,
            delta)
    else {   # delta needs to be formatted first
        ya <- array(y, c(M, J, T))
        yna <- apply(is.na(ya), c(1,3), all)
        delta <- formatDelta(delta, yna)
        out <- list(y=y, Xlam=Xlam, Xgam=Xgam, Xom=Xom, Xp=Xp, Xiota=Xiota,
                    Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
                    Xom.offset=Xom.offset, Xp.offset=Xp.offset,
                    Xiota.offset=Xiota.offset,
                    delta=delta, removed.sites=integer(0))
    }

    return(list(y = out$y, Xlam = out$Xlam, Xgam = out$Xgam,
                Xom = out$Xom, Xp = out$Xp, Xiota = out$Xiota,
                Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
                Xom.offset=Xom.offset, Xp.offset=Xp.offset,
                Xiota.offset=Xiota.offset, delta = out$delta,
                removed.sites = out$removed.sites, go.dims = go.dims))
})

#Need to do this hacky approach because class union of PCO and MMO doesn't work
#for reasons I don't understand
setMethod("getDesign", "unmarkedFrameMMO",
    function(umf, formula, na.rm=TRUE)
{
  class(umf)[1] <- "unmarkedFramePCO"
  getDesign(umf, formula, na.rm)
})

# need a getDesign for distsampOpen.... not sure how to set this up
# pcountOpenDS
setMethod("getDesign", "unmarkedFrameDSO",
    function(umf, formula, na.rm = TRUE)
{
    aschar1 <- as.character(formula)
    aschar2 <- as.character(formula[[2]])
    aschar3 <- as.character(formula[[2]][[2]])
    aschar4 <- as.character(formula[[2]][[2]][[2]])

    iotaformula <- as.formula(paste(aschar1[1], aschar1[3]))
    pformula <- as.formula(paste(aschar2[1], aschar2[3]))
    omformula <- as.formula(paste(aschar3[1], aschar3[3]))
    gamformula <- as.formula(paste(aschar4[1], aschar4[3]))
    lamformula <- as.formula(formula[[2]][[2]][[2]][[2]])

    y <- getY(umf)
    M <- nrow(y)
    T <- umf@numPrimary
    J <- ncol(y) / T
    delta <- umf@primaryPeriod

    if(is.null(umf@yearlySiteCovs))
        yearlySiteCovs <- data.frame(placeHolder = rep(1, M*T))
    else
        yearlySiteCovs <- umf@yearlySiteCovs

    ## add siteCovs in so they can be used as well
    if(!is.null(umf@siteCovs)) {
        sC <- umf@siteCovs[rep(1:M, each = T),,drop=FALSE]
        yearlySiteCovs <- cbind(yearlySiteCovs, sC)
        }

    if(is.null(siteCovs(umf)))
        siteCovs <- data.frame(placeHolder = rep(1, M))
    else
        siteCovs <- siteCovs(umf)

    Xlam.mf <- model.frame(lamformula, siteCovs, na.action = stats::na.pass)
    Xlam <- model.matrix(lamformula, Xlam.mf)
    Xlam.offset <- as.vector(model.offset(Xlam.mf))
    if(!is.null(Xlam.offset))
        Xlam.offset[is.na(Xlam.offset)] <- 0

    # Ignore last year of data
    transCovs <- yearlySiteCovs[-seq(T, M*T, by=T),,drop=FALSE]
    for(i in 1:ncol(transCovs))
        if(is.factor(transCovs[,i]))
            transCovs[,i] <- factor(transCovs[,i]) # drop unused levels

    Xiota.mf <- model.frame(iotaformula, transCovs, na.action = stats::na.pass)
    Xiota <- model.matrix(iotaformula, Xiota.mf)
    Xiota.offset <- as.vector(model.offset(Xiota.mf))
    if(!is.null(Xiota.offset))
        Xiota.offset[is.na(Xiota.offset)] <- 0

    #Detection uses yearlySiteCovs
    Xp.mf <- model.frame(pformula, yearlySiteCovs, na.action = stats::na.pass)
    Xp <- model.matrix(pformula, Xp.mf)
    Xp.offset <- as.vector(model.offset(Xp.mf))
    if(!is.null(Xp.offset))
        Xp.offset[is.na(Xp.offset)] <- 0

    Xgam.mf <- model.frame(gamformula, transCovs, na.action = stats::na.pass)
    Xgam <- model.matrix(gamformula, Xgam.mf)
    Xgam.offset <- as.vector(model.offset(Xgam.mf))
    if(!is.null(Xgam.offset))
        Xgam.offset[is.na(Xgam.offset)] <- 0

    Xom.mf <- model.frame(omformula, transCovs, na.action = stats::na.pass)
    Xom <- model.matrix(omformula, Xom.mf)
    Xom.offset <- as.vector(model.offset(Xom.mf))
    if(!is.null(Xom.offset))
        Xom.offset[is.na(Xom.offset)] <- 0

    if(is.null(Xlam.offset)) Xlam.offset <- rep(0, M)
    if(is.null(Xgam.offset)) Xgam.offset <- rep(0, M*(T-1))
    if(is.null(Xom.offset)) Xom.offset <- rep(0, M*(T-1))
    if(is.null(Xp.offset)) Xp.offset <- rep(0, M*T)
    if(is.null(Xiota.offset)) Xiota.offset<- rep(0, M*(T-1))

    # determine if gamma, omega, and iota are scalar, vector, or matrix valued
    # Runtime is much faster for scalars and vectors
    Xgo <- cbind(Xgam, Xom, Xiota)
    getGOdims <- function(x) {
        xm <- matrix(x, M, T-1, byrow=TRUE)
#        anyNA <- apply(is.na(xm), 1, any)
#        if(all(anyNA))
#            return("matrix")
#        xm <- xm[!anyNA,] # This is not 100% safe
#        nSites <- nrow(xm)
#        if(all(dim(unique(xm, MARGIN=1)) == c(1, T-1)))
#            return("rowvec")
#        else if(all(dim(unique(xm, MARGIN=2)) == c(nSites, 1)))
#            return("colvec")
#        else return("matrix")
        col.table <- apply(xm, 2, table)
        row.table <- apply(xm, 1, table)
        if(is.vector(col.table) & !is.list(col.table)) {
            return("rowvec")
        } else if(is.vector(row.table) & !is.list(row.table)) {
            return("colvec")
        } else
            return("matrix")
        }
    if(isTRUE(all.equal(gamformula,~1)) & isTRUE(all.equal(omformula, ~1)) &
      isTRUE(all.equal(iotaformula, ~1)))
        go.dims <- "scalar"
    else {
        go.dims.vec <- apply(Xgo, 2, getGOdims)
        if(all(go.dims.vec == "rowvec"))
            go.dims <- "rowvec"
        else if(all(go.dims.vec == "colvec"))
          ## NOTE: Temporary fix to the problem reported with
          ## time-only-varying covariates
            go.dims <- "matrix" ##"colvec"
        else
            go.dims <- "matrix"
    }

    if(na.rm){
      out <- handleNA(umf, Xlam, Xgam, Xom, Xp, Xiota,
                      Xlam.offset, Xgam.offset, Xom.offset, Xp.offset,
                      Xiota.offset, delta)
    } else {
      # delta needs to be formatted first
      ya <- array(y, c(M, J, T))
      yna <- apply(is.na(ya), c(1,3), all)
      delta <- formatDelta(delta, yna)
      out <- list(y=y, Xlam=Xlam, Xgam=Xgam, Xom=Xom, Xp=Xp, Xiota=Xiota,
                  Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
                  Xom.offset=Xom.offset, Xp.offset=Xp.offset,
                  Xiota.offset=Xiota.offset,
                  delta=delta, removed.sites=integer(0))
    }

    return(list(y = out$y, Xlam = out$Xlam, Xgam = out$Xgam,
                Xom = out$Xom, Xp = out$Xp, Xiota = out$Xiota,
                Xlam.offset=out$Xlam.offset, Xgam.offset=out$Xgam.offset,
                Xom.offset=out$Xom.offset, Xp.offset=out$Xp.offset,
                Xiota.offset=out$Xiota.offset, delta = out$delta,
                removed.sites = out$removed.sites, go.dims = go.dims))
})








#setMethod("handleNA", "unmarkedFramePCOorMMO",
setMethod("handleNA", "unmarkedFramePCO",
    function(umf, Xlam, Xgam, Xom, Xp, Xiota, Xlam.offset, Xgam.offset,
             Xom.offset, Xp.offset, Xiota.offset, delta)
{
	obsToY <- obsToY(umf)
	if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

	M <- numSites(umf)
	T <- umf@numPrimary
	y <- getY(umf)
	J <- ncol(y) / T
  R <- obsNum(umf)

	Xlam.long <- Xlam[rep(1:M, each = J*T),]
	Xlam.long.na <- is.na(Xlam.long)

	long.na <- function(x) {
            x.mat <- matrix(x, M, R, byrow = TRUE)
            x.mat <- is.na(x.mat)
            x.mat <- x.mat %*% obsToY
            x.long <- as.vector(t(x.mat))
            x.long > 0
        }

        o2y2 <- diag(T)
        o2y2 <- o2y2[-T, -T]

	long.na2 <- function(x) {
            x.mat <- matrix(x, M, T-1, byrow = TRUE)
            x.mat <- is.na(x.mat)
            x.mat <- x.mat %*% o2y2
            x.long <- as.vector(t(x.mat))
            x.long > 0
        }

	Xp.long.na <- apply(Xp, 2, long.na)
	Xp.long.na <- apply(Xp.long.na, 1, any)

        Xgam.long.na <- apply(Xgam, 2, long.na2)
	Xgam.long.na <- apply(Xgam.long.na, 1, any)
	Xom.long.na <- apply(Xom, 2, long.na2)
	Xom.long.na <- apply(Xom.long.na, 1, any)

	y.long <- as.vector(t(y))
	y.long.na <- is.na(y.long)

#  delta.long <- as.vector(t(delta))
#	delta.long.na <- is.na(delta.long)

	covs.na <- apply(cbind(Xlam.long.na, Xp.long.na), 1, any)
	covs.na2 <- apply(cbind(Xgam.long.na, Xom.long.na), 1, any)
	covs.na3 <- rep(covs.na2, each=J)
	# If gamma[1, 1] is NA, remove y[1, 2]
	#common <- 1:(M*J*(T-1))
        ignore <- rep(seq(1, M*J*T, by=J*T), each=J) + 0:(J-1)
        covs.na[-ignore] <- covs.na[-ignore] | covs.na3

	## are any NA in covs not in y already?
	y.new.na <- covs.na & !y.long.na

	if(sum(y.new.na) > 0) {
            y.long[y.new.na] <- NA
            warning("Some observations have been discarded because corresponding covariates were missing.", call. = FALSE)
        }

	y.wide <- matrix(y.long, nrow=M, ncol=J*T, byrow = TRUE)
#	delta <- matrix(delta.long, nrow=M, ncol=T, byrow = TRUE)
	sites.to.remove <- apply(y.wide, 1, function(x) all(is.na(x)))
  # Should also remove sites with no omega and gamma before an observation
  # remove all observations before the one after the first real omega/gamma
#  covs.na2.mat <- matrix(covs.na2, M, T-1, byrow=TRUE)
#  last.y <- apply(y, 1, function(x) max(which(!is.na(y))))
#  last.go <- apply(covs.na2.mat, 1, function(x)
#      if(any(x)) {
#          if(all(x))
#              return(0)
#          else
#              return(max(which(!x)))
#          }
#      else
#          return(T-1))
#  no.go.before.y <- last.y <= last.go

	ya <- array(y.wide, c(M, J, T))
	yna <- apply(is.na(ya), c(1,3), all)
	delta <- formatDelta(delta, yna)

	num.to.remove <- sum(sites.to.remove)
	if(num.to.remove > 0) {
            y.wide <- y.wide[!sites.to.remove, ,drop = FALSE]
            Xlam <- Xlam[!sites.to.remove, ,drop = FALSE]
            Xlam.offset <- Xlam.offset[!sites.to.remove]
            Xgam <- Xgam[!sites.to.remove[rep(1:M, each = T-1)],,
                         drop = FALSE]
            Xgam.offset <- Xgam.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                         drop = FALSE]
            Xom <- Xom[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            Xom.offset <- Xom.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            Xp <- Xp[!sites.to.remove[rep(1:M, each = J*T)],,
                     drop = FALSE]
            Xp.offset <- Xp.offset[!sites.to.remove[rep(1:M, each = J*T)],,
                     drop = FALSE]
            Xiota <- Xiota[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            Xiota.offset <- Xiota.offset[!sites.to.remove[rep(1:M, each = T-1)],,
                       drop = FALSE]
            delta <- delta[!sites.to.remove, ,drop =FALSE]
            warning(paste(num.to.remove, "sites have been discarded because of missing data."), call.=FALSE)
	}

	list(y = y.wide, Xlam = Xlam, Xgam = Xgam, Xom = Xom, Xp = Xp, Xiota = Xiota,
             Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
             Xom.offset=Xom.offset, Xp.offset=Xp.offset,
             Xiota.offset=Xiota.offset, delta = delta,
             removed.sites = which(sites.to.remove))
})


setMethod("handleNA", "unmarkedFrameDSO",
    function(umf, Xlam, Xgam, Xom, Xp, Xiota, Xlam.offset, Xgam.offset,
             Xom.offset, Xp.offset, Xiota.offset, delta)
{

  y <- getY(umf)
  M <- numSites(umf)
  T <- umf@numPrimary
  J <- ncol(y) / T

  ymat <- array(y, c(M,J,T))

  #Apply NAs for entire observation if any intervals are NA
  ymat <- array(y, c(M,J,T))
  obs_any_na <- apply(ymat, c(1,3), function(x) any(is.na(x)))
  any_na_ind <- which(obs_any_na, arr.ind=TRUE)
  if(sum(obs_any_na)>0){
    for(i in 1:nrow(any_na_ind)){
      ymat[any_na_ind[i,1], ,any_na_ind[i,2]] <- NA
    }
  }
  original_na_count <- sum(is.na(ymat))

  #NAs in lambda design matrix - have to delete entire site
  site_NA <- apply(Xlam, 1, function(x) any(is.na(x)))
  ymat[site_NA,,] <- NA

  find_na <- function(ymat, dm, M, T){
    dm <- matrix(apply(dm, 1, function(x) any(is.na(x))), M, T, byrow=TRUE)
    na_ind <- which(dm, arr.ind=TRUE)
    if(sum(dm)>0){
      for (i in 1:nrow(na_ind)){
        ymat[na_ind[i,1], ,na_ind[i,2]] <- NA
      }
    }
    ymat
  }

  #NAs in other design mats
  ymat <- find_na(ymat, Xgam, M, T-1)
  ymat <- find_na(ymat, Xom, M, T-1)
  ymat <- find_na(ymat, Xiota, M, T-1)
  ymat <- find_na(ymat, Xp, M, T)

  if(sum(is.na(ymat)) > original_na_count){
    warning("Some observations have been discarded because corresponding covariates were missing.", call. = FALSE)
  }

  #Format delta
  ytna <- apply(ymat, c(1,3), function(x) all(is.na(x)))
	delta <- formatDelta(delta, ytna)

  #Reformat y
  y <- matrix(ymat, nrow=M)

  #Remove sites
  sites.to.remove <- which(apply(ymat, 1, function(x) all(is.na(x))))
  nrem <- length(sites.to.remove)
  if(nrem > 0){
    y <- y[-sites.to.remove,]
    Xlam <- Xlam[-sites.to.remove,,drop=FALSE]
    Xlam.offset <- Xlam.offset[-sites.to.remove]
    delta <- delta[-sites.to.remove,,drop=FALSE]

    ysc_rem <- rep(1:M, each=(T-1)) %in% sites.to.remove
    Xgam <- Xgam[!ysc_rem,, drop=FALSE]
    Xgam.offset <- Xgam.offset[!ysc_rem]
    Xom <- Xom[!ysc_rem,, drop=FALSE]
    Xom.offset <- Xom.offset[!ysc_rem]
    Xiota <- Xiota[!ysc_rem,, drop=FALSE]
    Xiota.offset <- Xiota.offset[!ysc_rem]

    p_rem <- rep(1:M, each=T) %in% sites.to.remove
    Xp <- Xp[!p_rem,, drop=FALSE]
    Xp.offset <- Xp.offset[!p_rem]

    warning(paste(nrem, "sites have been discarded because of missing data."), call.=FALSE)
  }

	list(y = y, Xlam = Xlam, Xgam = Xgam, Xom = Xom, Xp = Xp, Xiota = Xiota,
             Xlam.offset=Xlam.offset, Xgam.offset=Xgam.offset,
             Xom.offset=Xom.offset, Xp.offset=Xp.offset,
             Xiota.offset=Xiota.offset, delta = delta,
             removed.sites = sites.to.remove)
})


# UnmarkedFrameGMN

setMethod("getDesign", "unmarkedFrameG3",
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
    R <- obsNum(umf) # 2*T for double observer sampling
                     # 1*T for distance sampling
                     # nPasses*T for removal sampling

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

    # add site and yearlysite covariates, which contain siteCovs
    cnames <- c(colnames(obsCovs), colnames(yearlySiteCovs))
    obsCovs <- cbind(obsCovs, yearlySiteCovs[rep(1:(M*T), each = R/T),])
    colnames(obsCovs) <- cnames

    # add observation number if not present
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

    return(list(y = out$y, Xlam = out$Xlam, Xphi = out$Xphi,
                Xdet = out$Xdet,
                Xlam.offset = out$Xlam.offset,
                Xphi.offset = out$Xphi.offset,
                Xdet.offset = out$Xdet.offset,
                removed.sites = out$removed.sites))
})





setMethod("handleNA", "unmarkedFrameG3",
    function(umf, Xlam, Xlam.offset, Xphi, Xphi.offset, Xdet, Xdet.offset)
{

#    browser()

    obsToY <- obsToY(umf)
    if(is.null(obsToY)) stop("obsToY cannot be NULL to clean data.")

    M <- numSites(umf)
    T <- umf@numPrimary
    R <- obsNum(umf)
    J <- numY(umf)/T

    # treat Xphi and Xlam together
    X <- cbind(Xphi, Xlam[rep(1:M, each = T), ])

    X.na <- is.na(X)
    X.long.na <- X.na[rep(1:(M*T), each = J),]

    Xdet.long.na <- apply(Xdet, 2, function(x) {
        x.mat <- matrix(x, M, R, byrow = TRUE)
        x.mat <- is.na(x.mat)
        x.mat <- x.mat %*% obsToY
        x.long <- as.vector(t(x.mat))
        x.long > 0
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
        y <- y[!sites.to.remove,, drop = FALSE]
        Xlam <- Xlam[!sites.to.remove,, drop = FALSE]
        Xlam.offset <- Xlam.offset[!sites.to.remove]
        Xphi <- Xphi[!sites.to.remove[rep(1:M, each = T)],, drop = FALSE]
        Xphi.offset <- Xphi.offset[!sites.to.remove[rep(1:M, each = T)]]
        Xdet <- Xdet[!sites.to.remove[rep(1:M, each = R)],,
                     drop=FALSE]
        Xdet.offset <- Xdet.offset[!sites.to.remove[rep(1:M, each=R)]]
        warning(paste(num.to.remove,
                      "sites have been discarded because of missing data."), call.=FALSE)
    }
    list(y = y, Xlam = Xlam, Xlam.offset = Xlam.offset, Xphi = Xphi,
        Xphi.offset = Xphi.offset, Xdet = Xdet, Xdet.offset = Xdet.offset,
        removed.sites = which(sites.to.remove))
})

