# Fit the occupancy model COP
# (Counting Occurrences Process)

# Occupancy
# z_i ~ Bernoulli(psi_i)
# 
# with:
#   z_i = Occupancy state of site i
#       = 1 if the site i is occupied
#       = 0 else
#   psi_i = Occupancy probability of site i

# Detection
# N_ij | z_i = 1 ~ Poisson(lambda_ij*L_ij)
# N_ij | z_i = 0 ~ 0
# 
# with:
#   N_ij = Number of detections of site i during observation j
#   z_i = Occupancy state of site i
#   lambda_ij = Detection rate of the observation j in site i
#   L_ij = Length/Duration of the observation j in site i

# CLASSES ----------------------------------------------------------------------

## unmarkedFrameOccuCOP class ----
setClass(
  "unmarkedFrameOccuCOP",
  representation(L = "matrix"),
  contains = "unmarkedFrame",
  validity = function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    y_integers = sapply(object@y, check.integer, na.ignore = T)
    if (!all(y_integers)) {
      errors <- c(errors,
                  paste(
                    "Count detection should be integers. Non-integer values:",
                    paste(object@y[which(!y_integers)], collapse = ', ')
                  ))
    }
    if (!all(all(dim(object@L) == dim(object@y)))){
      errors <- c( errors, paste(
        "L should be a matrix of the same dimension as y, with M =", M,
        "rows (sites) and J =", J, "columns (sampling occasions)."
      ))}
    if (length(errors) == 0) TRUE
    else errors
  }
)

## unmarkedFitOccuCOP class ----
setClass("unmarkedFitOccuCOP",
         representation(removed_obs = "matrix",
                        formlist = "list"),
         contains = "unmarkedFit")


# Methods ----------------------------------------------------------------------

## getDesign method ----
setMethod(
  "getDesign", "unmarkedFrameOccuCOP",
  function(umf, formlist, na.rm = TRUE) {
    
    "
    getDesign convert the information in the unmarkedFrame to a format
    usable by the likelihood function:
    - It creates model design matrices for fixed effects (X) for each parameter (here lambda, psi) 
    - It creates model design matrices for random effects (Z) for each parameter (here lambda, psi)
    - It handles missing data
    "
    
    # Retrieve useful informations from umf
    M <- numSites(umf)
    J <- obsNum(umf)
    y <- getY(umf)
    L <- getL(umf)
    
    # Occupancy submodel -------------------------------------------------------
    # Retrieve the fixed-effects part of the formula
    psiformula <- remove_bars(as.formula(formlist$psiformula))
    psiVars <- all.vars(psiformula)
    
    # Retrieve the site covariates
    sc <- siteCovs(umf)
    if(is.null(sc)) {
      sc <- data.frame(.dummy = rep(0, M))
    }
    
    # Check for missing variables
    psiMissingVars <- psiVars[!(psiVars %in% names(sc))]
    if (length(psiMissingVars) > 0) {
      stop(paste0(
        "Variable(s) '",
        paste(psiMissingVars, collapse = "', '"),
        "' not found in siteCovs"
      ), call. = FALSE)
    }

    # State model matrix for fixed effects
    Xpsi <- model.matrix(
      psiformula,
      model.frame(psiformula, sc, na.action = NULL)
    )
    # State model matrix for random effects
    Zpsi <- get_Z(formlist$psiformula, sc)
    
    # Detection submodel -------------------------------------------------------
    
    # Retrieve the fixed-effects part of the formula
    lambdaformula <- remove_bars(as.formula(formlist$lambdaformula))
    lambdaVars <- all.vars(lambdaformula)
    
    # Retrieve the observation covariates
    oc <- obsCovs(umf)
    if(is.null(oc)) {
      oc <- data.frame(.dummy = rep(0, M*J))
    }
    
    # Check for missing variables
    lambdaMissingVars <- lambdaVars[!(lambdaVars %in% names(oc))]
    if (length(lambdaMissingVars) > 0) {
      stop(paste(
        "Variable(s)",
        paste(lambdaMissingVars, collapse = ", "),
        "not found in obsCovs"
      ), call. = FALSE)
    }
    
    # Detection model matrix for fixed effects
    Xlambda <- model.matrix(
      lambdaformula,
      model.frame(lambdaformula, oc, na.action = NULL)
    )
    # Detection model matrix for random effects
    Zlambda <- get_Z(formlist$lambdaformula, oc)
    
    # Missing data -------------------------------------------------------------
    
    # Missing count data
    missing_y <- is.na(y)
    
    # Missing site covariates
    # (TRUE if at least one site covariate is missing in a site)
    missing_sc <- apply(Xpsi, 1, function(x) any(is.na(x)))
    
    # Missing observation covariates
    # (TRUE if at least one observation covariate is missing in a sampling occasion in a site)
    missing_oc <- apply(Xlambda, 1, function(x) any(is.na(x)))
    
    # Matrix MxJ of values to not use because there is some data missing
    removed_obs <- 
      # If there is count data missing in site i and obs j
      missing_y | 
      # If there is any site covariate missing in site i
      replicate(n = J, missing_sc) |
      # If there is any observation covariate missing in site i and obs j
      matrix(missing_oc, M, J, byrow = T)
    
    if (any(removed_obs)) {
      if (na.rm) {
        nb_missing_sites <- sum(rowSums(!removed_obs) == 0)
        nb_missing_observations <- sum(is.na(removed_obs))
        warning("There is missing data: ",
                sum(missing_y), " missing count data, ",
                sum(missing_sc), " missing site covariate(s), ",
                sum(missing_oc), " missing observation covariate(s). ",
                "Data from only ", (M*J)-sum(removed_obs), " observations out of ", (M*J), " are used, ",
                "from ", M-nb_missing_sites, " sites out of ", M, ".\n\t"
        )
      } else {
        stop("na.rm=FALSE and there is missing data :\n\t",
             sum(missing_y), " missing count data (y)\n\t",
             sum(missing_sc), " missing site covariates (siteCovs)\n\t",
             sum(missing_oc), " missing observation covariates (obsCovs)")
      }
    }
    
    # Output -------------------------------------------------------------------
    return(list(
      y = y,
      Xpsi = Xpsi,
      Zpsi = Zpsi,
      Xlambda = Xlambda,
      Zlambda = Zlambda,
      removed_obs = removed_obs
    ))
  })


## getL method ----
setGeneric("getL", function(object) standardGeneric("getL"))
setMethod("getL", "unmarkedFrameOccuCOP", function(object) {
  return(object@L)
})


## show method ----
setMethod("show", "unmarkedFrameOccuCOP", function(object) {
  J <- ncol(object@L)
  df_unmarkedFrame <- as(object, "data.frame")
  df_L <- data.frame(object@L)
  colnames(df_L) <- paste0("L.", 1:J)
  if (ncol(df_unmarkedFrame) > J) {
    df <- cbind(df_unmarkedFrame[, 1:J, drop = FALSE], 
                df_L, 
                df_unmarkedFrame[, (J + 1):ncol(df_unmarkedFrame), drop = FALSE])
  } else {
    df <- cbind(df_unmarkedFrame[, 1:J], 
                df_L)
  }
  cat("Data frame representation of unmarkedFrame object.\n")
  print(df)
})
# LP note: as is defined in unmarkedFrame.R part "COERCION"


## summary method ----
setMethod("summary", "unmarkedFrameOccuCOP", function(object,...) {
  cat("unmarkedFrameOccuCOP Object\n\n")
  
  cat(nrow(object@y), "sites\n")
  cat("Maximum number of sampling occasions per site:",obsNum(object),"\n")
  mean.obs <- mean(rowSums(!is.na(getY(object))))
  cat("Mean number of sampling occasions per site:",round(mean.obs,2),"\n")
  cat("Sites with at least one detection:",
      sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))),
      "\n\n")
  
  cat("Tabulation of y observations:")
  print(table(object@y, exclude=NULL))
  
  if(!is.null(object@L)) {
    cat("\nTabulation of sampling occasions length:")
    print(table(object@L))
  }

  if(!is.null(object@siteCovs)) {
    cat("\nSite-level covariates:\n")
    print(summary(object@siteCovs))
  }
  
  if(!is.null(object@obsCovs)) {
    cat("\nObservation-level covariates:\n")
    print(summary(object@obsCovs))
  }
})


## umf[i, j] ----
setMethod("[", c("unmarkedFrameOccuCOP", "numeric", "numeric", "missing"),
          function(x, i, j) {
            # Gey dimensions of x
            M <- numSites(x)
            J <- obsNum(x)
            
            if (length(i) == 0 & length(j) == 0) {
              return(x)
            }
            
            # Check i
            if (any(i < 0) &&
                any(i > 0)) {
              stop("i must be all positive or all negative indices.")
            }
            if (all(i < 0)) {
              i <- (1:M)[i]
            }
            
            # Check j
            if (any(j < 0) &&
                any(j > 0)) {
              stop("j must be all positive or all negative indices.")
            }
            if (all(j < 0)) {
              j <- (1:J)[j]
            }
            
            # y observation count data subset
            y <- getY(x)[i, j, drop = FALSE]
            if (min(length(i), length(j)) == 1) {
              y <- t(y)
            }
            
            # L subset
            L <- x@L[i, j, drop = FALSE]
            if (min(length(i), length(j)) == 1) {
              L <- t(L)
            }
            
            # siteCovs subset
            siteCovs <- siteCovs(x)
            if (!is.null(siteCovs)) {
              siteCovs <- siteCovs(x)[i, , drop = FALSE]
            }
            
            # obsCovs subset
            obsCovs <- obsCovs(x)
            if (!is.null(obsCovs)) {
              MJ_site <- rep(1:M, each = J)
              MJ_obs <- rep(1:J, times = M)
              obsCovs <- obsCovs[((MJ_obs %in% j) & (MJ_site %in% i)), , drop = FALSE]
              rownames(obsCovs) <- NULL
            }
            
            # Recreate umf
            new(
              Class = "unmarkedFrameOccuCOP",
              y = y,
              L = L,
              siteCovs = siteCovs,
              obsCovs = obsCovs,
              obsToY = diag(length(j)),
              mapInfo = x@mapInfo
            )
          })


## umf[i, ] ----
setMethod("[", c("unmarkedFrameOccuCOP", "numeric", "missing", "missing"),
          function(x, i) {
            x[i, 1:obsNum(x)]
          })

## umf[, j] ----
setMethod("[", c("unmarkedFrameOccuCOP", "missing", "numeric", "missing"),
          function(x, j) {
            x[1:numSites(x), j]
          })


## fl_getY ----
setMethod("fl_getY", "unmarkedFitOccuCOP", function(fit, ...){
  getDesign(getData(fit), fit@formlist)$y
})


## predict_inputs_from_umf ----
setMethod("predict_inputs_from_umf", "unmarkedFitOccuCOP",
          function(object, type, newdata, na.rm, re.form = NULL) {
            designMats = getDesign(umf = newdata,
                                   formlist = object@formlist,
                                   na.rm = na.rm)
            if (type == "psi") list_els <- c("Xpsi", "Zpsi")
            if (type == "lambda") list_els <- c("Xlambda", "Zlambda")
            X <- designMats[[list_els[1]]]
            if (is.null(re.form)) X <- cbind(X, designMats[[list_els[2]]])
            return(list(X = X, offset = NULL))
          })


## get_formula ----
setMethod("get_formula", "unmarkedFitOccuCOP", function(object, type, ...) {
  fl <- object@formlist
  switch(type, psi = fl$psiformula, lambda = fl$lambdaformula)
})


## get_orig_data ----
setMethod("get_orig_data", "unmarkedFitOccuCOP", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=FALSE)
  datatype <- switch(type, psi = 'site_covs', lambda = 'obs_covs')
  clean_covs[[datatype]]
})


## getP ----
setMethod("getP", "unmarkedFitOccuCOP", function(object, na.rm = TRUE) {
  data <- object@data
  M = nrow(getY(data))
  J = ncol(getY(data))
  des <- getDesign(data, object@formlist, na.rm = na.rm)
  matLambda =  do.call(object["lambda"]@invlink, 
                       list(matrix(
                         as.numeric(des$Xlambda %*% coef(object, 'lambda')),
                         nrow = M, ncol = J, byrow = T)))
  return(matLambda)
})


## fitted ----
setMethod("fitted", "unmarkedFitOccuCOP", function(object, na.rm = FALSE) {
  data <- object@data
  M = nrow(getY(data))
  J = ncol(getY(data))
  des <- getDesign(data, object@formlist, na.rm = na.rm)
  estim_psi = as.numeric(do.call(object["psi"]@invlink,
                                 list(as.matrix(des$Xpsi %*% coef(object, 'psi')))))
  estim_lambda = do.call(object["lambda"]@invlink, 
                         list(matrix(
                           as.numeric(des$Xlambda %*% coef(object, 'lambda')),
                           nrow = M, ncol = J, byrow = T)))
  return(estim_psi * estim_lambda)
})


## residuals ----
setMethod("residuals", "unmarkedFitOccuCOP", function(object) {
  y <- getY(object@data)
  e <- fitted(object)
  r <- y - e
  return(r)
})


## plot ----
setMethod("plot", c(x = "unmarkedFitOccuCOP", y = "missing"), function(x, y, ...) {
  y <- getY(x)
  r <- residuals(x)
  e <- fitted(x)
  
  old_graph <- graphics::par("mfrow", "mar")
  on.exit(graphics::par(mfrow = old_graph$mfrow, mar = old_graph$mar))
  
  {
    graphics::par(mfrow = c(2, 1))
    graphics::par(mar = c(0, 5, 2, 2))
    plot(e, y,
         ylab = "Observed data",
         xlab = "Predicted data",
         xaxt = 'n')
    abline(a = 0, b = 1, lty = 3, col = "red")
    title(main = "COP model - detection events count", outer = F)
    
    graphics::par(mar = c(4, 5, 0.5, 2))
    plot(e, r,
         ylab = "Residuals",
         xlab = "Predicted data")
    abline(h = 0, lty = 3, col = "red")
  }
})


## get_umf_components ----
setMethod("get_umf_components", "unmarkedFitOccuCOP",
  function(object, formulas, guide, design, ...){
    sc <- generate_data(formulas$psi, guide, design$M)
    oc <- generate_data(formulas$lambda, guide, design$J*design$M)
    yblank <- matrix(0, design$M, design$J)
    list(y=yblank, siteCovs=sc, obsCovs=oc)
})


## simulate_fit ----
setMethod("simulate_fit", "unmarkedFitOccuCOP",
  function(object, formulas, guide, design, ...){
    # Generate covariates and create a y matrix of zeros
    parts <- get_umf_components(object, formulas, guide, design, ...)
    umf <- unmarkedFrameOccuCOP(y = parts$y, siteCovs = parts$siteCovs, obsCovs=parts$obsCovs)
    fit <- suppressMessages(
      occuCOP(
        data = umf,
        psiformula = formula(formulas$psi),
        lambdaformula = formula(formulas$lambda),
        se = FALSE,
        control = list(maxit = 1)
      )
    )
    return(fit)
})


## simulate ----
setMethod("simulate", "unmarkedFitOccuCOP",
  function(object, nsim = 1, seed = NULL, na.rm = TRUE){
  # set.seed(seed) 
  # Purposefully not implemented
  formula <- object@formula
  umf <- object@data
  designMats <- getDesign(umf = umf, formlist = object@formlist, na.rm = na.rm)
  y <- designMats$y
  M <- nrow(y)
  J <- ncol(y)
  
  # Occupancy probability psi depending on the site covariates
  psiParms = coef(object, type = "psi", fixedOnly = FALSE)
  psi <- as.numeric(plogis(as.matrix(designMats$Xpsi %*% psiParms)))
  
  # Detection rate lambda depending on the observation covariates
  lambda = getP(object = object)
  
  # Simulations
  simList <- vector("list", nsim)
  for(i in 1:nsim) {
    Z <- rbinom(M, 1, psi)
    # Z[object@knownOcc] <- 1
    y = matrix(rpois(n = M * J, lambda = as.numeric(t(lambda))),
               nrow = M, ncol = J, byrow = T) * Z
    simList[[i]] <- y
  }
  return(simList)
})


## nonparboot ----
setMethod("nonparboot", "unmarkedFitOccuCOP",
  function(object, B = 0, keepOldSamples = TRUE, ...) {
  stop("Not currently supported for unmarkedFitOccuCOP", call.=FALSE)
})


## ranef ----
setMethod("ranef", "unmarkedFitOccuCOP", function(object, ...) {
  # Sites removed (srm) and sites kept (sk)
  srm <- object@sitesRemoved
  if (length(srm) > 0) {
    sk = 1:numSites(getData(object))[-srm]
  } else{
    sk = 1:numSites(getData(object))
  }
  
  # unmarkedFrame informations
  M <- length(sk)
  J <- obsNum(getData(object))
  y <- getY(getData(object))[sk,]
  
  # Estimated parameters
  psi <- predict(object, type = "psi")[sk, 1]
  lambda <- getP(object)[sk,]
  
  # Estimate posterior distributions
  z = c(0, 1)
  post <- array(0, c(M, 2, 1), dimnames = list(NULL, z))
  for (i in 1:M) {
    # psi density
    f <- dbinom(x = z,
                size = 1,
                prob = psi[i])
    
    # lambda density
    g <- c(1, 1)
    for (j in 1:J) {
      if (is.na(y[i, j]) | is.na(lambda[i, j])) {
        next
      }
      g = g * dpois(x = y[i, j], lambda = lambda[i, j] * z)
    }
    
    # psi*lambda density
    fudge <- f * g
    post[i, , 1] <- fudge / sum(fudge)
  }
  
  new("unmarkedRanef", post = post)
})


# Useful functions -------------------------------------------------------------

check.integer <- function(x, na.ignore = F) {
  if (is.na(x) & na.ignore) {
    return(T)
  }
  x %% 1 == 0
}

# unmarkedFrame ----------------------------------------------------------------

unmarkedFrameOccuCOP <- function(y, L, siteCovs = NULL, obsCovs = NULL) {
  
  # Verification that they are non-NA data in y
  if (all(is.na(y))) {
    stop("y only contains NA. y needs to contain non-NA integers.")
  }
  
  # Verification that these are count data (and not detection/non-detection)
  if (max(y, na.rm = T) == 1) {
    warning("unmarkedFrameOccuCOP is for count data. ",
            "The data furnished appear to be detection/non-detection.")
  }
  
  # Number of sampling occasions
  J <- ncol(y)
  
  # If missing L: replace by matrix of 1
  # and the lambda will be the detection rate per observation length
  if (missing(L)) {
    L <- y * 0 + 1
    warning("L is missing, replacing it by a matrix of 1.")
  } else if (is.null(L)) {
    L <- y * 0 + 1
    warning("L is missing, replacing it by a matrix of 1.")
  }
  
  # Transformation observation covariates
  obsCovs <- covsToDF(
    covs = obsCovs,
    name = "obsCovs",
    obsNum = J,
    numSites = nrow(y)
  )
  
  # Create S4 object of class unmarkedFrameOccuCOP
  umf <- new(
    Class = "unmarkedFrameOccuCOP",
    y = y,
    L = L,
    siteCovs = siteCovs,
    obsCovs = obsCovs,
    obsToY = diag(J)
  )
  
  return(umf)
}


# occuCOP ----------------------------------------------------------------------

occuCOP <- function(data,
                    psiformula = ~1,
                    lambdaformula = ~1,
                    psistarts,
                    lambdastarts,
                    starts,
                    method = "BFGS",
                    se = TRUE,
                    engine = c("C", "R"),
                    na.rm = TRUE,
                    return.negloglik = NULL,
                    L1 = FALSE,
                    ...) {
  #TODO: random effects
  
  # Neg loglikelihood COP ------------------------------------------------------
  R_nll_occuCOP <- function(params) {
    
    # Reading and transforming params
    # Taking into account the covariates
    
    # psi as a function of covariates
    #   psi in params are transformed with a logit transformation (qlogis)
    #       so they're back-transformed to a proba with inverse logit (plogis)
    #   Xpsi is the matrix with occupancy covariates
    #   params is the vector with all the parameters
    #   psiIdx is the index of Occupancy Parameters in params
    psi <- plogis(Xpsi %*% params[psiIdx])
    
    # lambda as a function of covariates
    #   lambda in params are transformed with a log-transformation
    #          so they're back-transformed to a rate with exp here
    #   Xlambda is the matrix with detection covariates
    #   params is the vector with all the parameters
    #   lambdaIdx is the index of Occupancy Parameters in params
    lambda <- exp(Xlambda %*% params[lambdaIdx])
    
    # Listing sites analysed = sites not removed (due to NAs)
    if (length(sitesRemoved) > 0) {
      siteAnalysed = (1:M)[-sitesRemoved]
    } else {
      siteAnalysed = (1:M)
    }
    
    # Probability for each site (i)
    iProb <- rep(NA, M)
    
    for (i in siteAnalysed) {
      # iIdx is the index to access the vectorised vectors of all obs in site i
      iIdxall <- ((i - 1) * J + 1):(i * J)
      
      # Removing NAs
      iIdx = iIdxall[!removed_obsvec[iIdxall]]
      
      if (SitesWithDetec[i]) {
        # If there is at least one detection in site i
        iProb[i] = psi[i] * (
          (sum(lambda[iIdx] * Lvec[iIdx])) ^ sum(yvec[iIdx]) / 
            factorial(sum(yvec[iIdx])) *
            exp(-sum(lambda[iIdx] * Lvec[iIdx]))
        )
        
      } else {
        # If there is zero detection in site i
        iProb[i] = psi[i] * exp(-sum(lambda[iIdx] * Lvec[iIdx])) + (1 - psi[i]) 
        
      }
    }
    
    # log-likelihood
    ll = sum(log(iProb[siteAnalysed]))
    return(-ll)
  }
  
  # Check arguments ------------------------------------------------------------
  if (!is(data, "unmarkedFrameOccuCOP")) {
    stop("Data is not an unmarkedFrameOccuCOP object. See ?unmarkedFrameOccuCOP if necessary.")
  }
  stopifnot(class(psiformula) == "formula")
  stopifnot(class(lambdaformula) == "formula")
  if(!missing(psistarts)){stopifnot(class(psistarts) %in%  c("numeric", "double", "integer"))}
  if(!missing(lambdastarts)){stopifnot(class(lambdastarts) %in%  c("numeric", "double", "integer"))}
  se = as.logical(match.arg(
    arg = as.character(se),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  na.rm = as.logical(match.arg(
    arg = as.character(na.rm),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  engine <- match.arg(engine, c("C", "R"))
  L1 = as.logical(match.arg(
    arg = as.character(L1),
    choices = c("TRUE", "FALSE", "0", "1")
  ))
  
  
  # Do not yet manage random effects!!!
  if (has_random(psiformula) | has_random(lambdaformula)) {
    stop("occuCOP does not currently handle random effects.")
  }
  
  # Format input data ----------------------------------------------------------
  
  # Retrieve formlist
  formlist <- mget(c("psiformula", "lambdaformula"))
  
  # Get the design matrix (calling the getDesign method for unmarkedFrame)
  # For more informations, see: getMethod("getDesign", "unmarkedFrameOccuCOP")
  designMats <- getDesign(umf = data, formlist = formlist, na.rm = na.rm)
  
  # y is the count detection data (matrix of size M sites x J observations)
  y <- getY(data)
  
  # L is the length of observations (matrix of size M sites x J observations)
  L <- getL(data)
  if (!L1) {
    if (!any(is.na(L))) {
      if (all(L == 1)) {
        warning(
          "All observations lengths (L) are set to 1. ",
          "If they were not user-defined, lambda corresponds to the ",
          "detection rate multiplied by the observation length, ",
          "not just the detection rate per time-unit or space-unit.\n",
          "You can remove this warning by adding 'L1=TRUE' in the function inputs."
        )
      }
    }
  }
  
  # Xpsi is the fixed effects design matrix for occupancy
  Xpsi <- designMats$Xpsi
  
  # Xlambda is the fixed effects design matrix for detection rate
  Xlambda <- designMats$Xlambda
  
  # Zpsi is the random effects design matrix for occupancy
  Zpsi <- designMats$Zpsi
  
  # Zlambda is the random effects design matrix for detection rate
  Zlambda <- designMats$Zlambda
  
  # removed_obs is a M x J matrix of the observations removed from the analysis
  removed_obs <- designMats$removed_obs
  sitesRemoved <- unname(which(apply(removed_obs, 1, function(x) all(x))))
  
  # Number of sites
  M <- nrow(y)
  
  # Number of sampling occasions
  J <- ncol(y)
  
  # Total number of detection per site
  NbDetecPerSite = rowSums(y, na.rm=T)
  
  # Sites where there was at least one detection
  SitesWithDetec = NbDetecPerSite > 0
  
  # Number of sites where there was at least one detection
  NbSitesWithDetec = sum(SitesWithDetec)
  
  # Set up parameter names and indices-----------------------------------------
  
  # ParamPsi Occupancy parameter names
  ParamPsi <- colnames(Xpsi)
  
  # ParamLambda Detection parameter names
  ParamLambda <- colnames(Xlambda)
  
  # NbParamPsi Number of occupancy parameters
  NbParamPsi <- ncol(Xpsi)
  
  # NbParamLambda Number of detection parameters
  NbParamLambda <- ncol(Xlambda)
  
  # nP Total number of parameters
  nP <- NbParamPsi + NbParamLambda
  
  # psiIdx Index of the occupancy parameters
  psiIdx <- 1:NbParamPsi
  
  # lambdaIdx Index of the detection parameters
  lambdaIdx <- (NbParamPsi+1):nP
  
  # Re-format some variables for C and R engines
  #   From Matrix of dim MxJ to vector of length MxJ:
  #   c(ySite1Obs1, ySite1Obs2, ..., ySite1ObsJ, ysite2Obs1, ...)
  yvec <- as.numeric(t(y))
  Lvec <- as.numeric(t(L))
  removed_obsvec <- as.logical(t(removed_obs))
  
  # return.negloglik -----------------------------------------------------------
  if (!is.null(return.negloglik)) {
    df_NLL = data.frame(t(as.data.frame(return.negloglik)))
    rownames(df_NLL) = NULL
    colnames(df_NLL) = c(paste0("logit(psi).", ParamPsi),
                         paste0("log(lambda).", ParamLambda))
    df_NLL$nll = NA
    for (i in 1:nrow(df_NLL)) {
      df_NLL$nll[i] = R_nll_occuCOP(params = as.numeric(as.vector(df_NLL[i, -ncol(df_NLL)])))
    }
    return(df_NLL)
  }
  
  # nll function depending on engine -------------------------------------------
  if (identical(engine, "C")) {
    nll <- function(params) {
      nll_occuCOP(
        y = yvec,
        L = Lvec,
        Xpsi = Xpsi,
        Xlambda = Xlambda,
        beta_psi = params[psiIdx],
        beta_lambda = params[lambdaIdx],
        removed = removed_obsvec
      )
    }
  } else if (identical(engine, "R")) {
    nll <- R_nll_occuCOP
  }
  
  
  # Optimisation ---------------------------------------------------------------
  
  ## Checking the starting point for optim ----
  # Check if either (psistarts AND lambdastarts) OR starts is provided
  if (!missing(psistarts) & !missing(lambdastarts)) {
    # Both psistarts and lambdastarts provided
    if (!missing(starts)){
      if (!all(c(psistarts, lambdastarts) == starts)) {
        warning(
          "You provided psistarts, lambdastarts and starts. ",
          "Please provide either (psistarts AND lambdastarts) OR starts. ",
          "Using psistarts and lambdastarts."
        )
      }
    }
    if (length(lambdastarts) != NbParamLambda) {
      stop("lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
           "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula)
    }
    if (length(psistarts) != NbParamPsi) {
      stop("psistarts (", paste(psistarts, collapse = ", "), ") ",
           "should be of length ", NbParamPsi, " with psiformula ", psiformula)
    }
    starts <- c(psistarts, lambdastarts)
  } else if (!missing(starts)) {
    # starts provided
    if (length(starts) != nP) {
      stop("starts (", paste(starts, collapse = ", "), ") ",
           "should be of length ", nP, 
           " with psiformula ", psiformula, 
           " and lambdaformula ", lambdaformula)
    }
    
    psistarts <- starts[1:NbParamPsi]
    lambdastarts <- starts[(NbParamPsi + 1):(NbParamPsi + NbParamLambda)]
    
  } else {
    # No arguments provided, apply default values
    
    if (missing(lambdastarts)) {
      # If lambda starts argument was not given:
      # 0 by default 
      # so lambda = exp(0) = 1 by default
      lambdastarts = rep(0, NbParamLambda)
    } else if (length(lambdastarts) != NbParamLambda) {
      stop("lambdastarts (", paste(lambdastarts, collapse = ", "), ") ",
           "should be of length ", NbParamLambda, " with lambdaformula ", lambdaformula)
    }
    
    if (missing(psistarts)) {
      # If psi starts argument was not given
      # 0 by default 
      # so psi = plogis(0) = 0.5 by default
      psistarts = rep(0, NbParamPsi)
    } else if (length(psistarts) != NbParamPsi) {
      stop("psistarts (", paste(psistarts, collapse = ", "), ") ",
           "should be of length ", NbParamPsi, " with psiformula ", psiformula)
    }
    
    starts <- c(psistarts, lambdastarts)
  }
  
  ## Run optim ----
  opt <- optim(
    starts,
    nll,
    method = method,
    hessian = se,
    ...
  )
  
  # Get output -----------------------------------------------------------------
  covMat <- invertHessian(opt, nP, se)
  ests <- opt$par
  tmb_mod <- NULL
  fmAIC <- 2 * opt$value + 2 * nP
  
  # Organize effect estimates
  names(ests) <- c(ParamPsi, ParamLambda)
  psi_coef <- list(ests = ests[psiIdx], cov = as.matrix(covMat[psiIdx, psiIdx]))
  lambda_coef <- list(ests = ests[lambdaIdx],
                      cov = as.matrix(covMat[lambdaIdx, lambdaIdx]))
  
  # No random effects
  psi_rand_info <- lambda_rand_info <- list()
  
  # Create unmarkedEstimates ---------------------------------------------------
  psi_uE <- unmarkedEstimate(
    name = "Occupancy probability",
    short.name = "psi",
    estimates = psi_coef$ests,
    covMat = psi_coef$cov,
    fixed = 1:NbParamPsi,
    invlink = "logistic",
    invlinkGrad = "logistic.grad",
    randomVarInfo = psi_rand_info
  )
  
  lambda_uE <- unmarkedEstimate(
    name = "Detection rate",
    short.name = "lambda",
    estimates = lambda_coef$ests,
    covMat = lambda_coef$cov,
    fixed = 1:NbParamLambda,
    invlink = "exp",
    invlinkGrad = "exp",
    randomVarInfo = lambda_rand_info
  )
  
  estimateList <- unmarkedEstimateList(list(psi = psi_uE, lambda = lambda_uE))
  
  # Create unmarkedFit object--------------------------------------------------
  umfit <- new(
    "unmarkedFitOccuCOP",
    fitType = "occuCOP",
    call = match.call(),
    formula = as.formula(paste(
      formlist["lambdaformula"], formlist["psiformula"], collapse = ""
    )),
    formlist = formlist,
    data = data,
    estimates = estimateList,
    sitesRemoved = sitesRemoved,
    removed_obs = removed_obs,
    AIC = fmAIC,
    opt = opt,
    negLogLike = opt$value,
    nllFun = nll,
    TMB = tmb_mod
  )
  
  return(umfit)
}
