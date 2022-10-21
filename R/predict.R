# General predict method-------------------------------------------------------

# Common predict function for all fit types
# with exception of occuMulti and occuMS (at the end of this file)
setMethod("predict", "unmarkedFit",
  function(object, type, newdata, backTransform = TRUE, na.rm = TRUE,
           appendData = FALSE, level=0.95, re.form=NULL, ...){

  # If no newdata, get actual data
  if(missing(newdata) || is.null(newdata)) newdata <- object@data

  # Check inputs
  check_predict_arguments(object, type, newdata)

  # Get model matrix (X) and offset
  # If newdata is an unmarkedFrame, use getDesign via predict_inputs_from_umf()
  is_raster <- FALSE
  if(inherits(newdata, "unmarkedFrame")){
    # Generate model matrix and offsets
    pred_inps <- predict_inputs_from_umf(object, type, newdata, na.rm, re.form)
  } else {
    # If newdata is provided
    # 1. Get original data and appropriate formula for type
    orig_data <- get_orig_data(object, type)
    orig_formula <- get_formula(object, type)

    # 2. If newdata is raster, get newdata from raster as data.frame
    if(inherits(newdata, c("RasterLayer","RasterStack"))){
      if(!require(raster)) stop("raster package required", call.=FALSE)
      is_raster <- TRUE
      orig_raster <- newdata
      newdata <- newdata_from_raster(newdata, all.vars(orig_formula))
    }

    # 3. Make model matrix and offset with newdata, informed by original data
    pred_inps <- make_mod_matrix(orig_formula, orig_data, newdata, re.form)
  }

  # Calculate predicted values in chunks (for speed) based on X and offset
  out <- predict_by_chunk(object, type, level, pred_inps$X, pred_inps$offset,
                          chunk_size = 70, backTransform, re.form)

  # Convert output to raster if newdata was raster
  if(is_raster){
    out <- raster_from_predict(out, orig_raster, appendData)
  } else if(appendData){
    # Append data if needed
    out <- data.frame(out, as(newdata, "data.frame"))
  }

  out

})

# Function to make model matrix and offset from formula, newdata and original data
# This function makes sure factor levels in newdata match, and that
# any functions in the formula are handled properly (e.g. scale)
make_mod_matrix <- function(formula, data, newdata, re.form=NULL){
  form_nobars <- lme4::nobars(formula)
  mf <- model.frame(form_nobars, data, na.action=stats::na.pass)
  X.terms <- stats::terms(mf)
  fac_cols <- data[, sapply(data, is.factor), drop=FALSE]
  xlevs <- lapply(fac_cols, levels)
  xlevs <- xlevs[names(xlevs) %in% names(mf)]
  nmf <- model.frame(X.terms, newdata, na.action=stats::na.pass, xlev=xlevs)
  #X <- model.matrix(X.terms, newdata, xlev=xlevs)
  X <- model.matrix(form_nobars, nmf)
  offset <- model.offset(nmf)
  if(is.null(re.form) & !is.null(lme4::findbars(formula))){
    Z <- get_Z(formula, data, newdata)
    X <- cbind(X, Z)
  }
  list(X=X, offset=offset)
}


# Fit-type specific methods----------------------------------------------------

# Fit type-specific methods to generate different components of prediction
# 1. check_predict_arguments(): Check arguments
# 2. predict_inputs_from_umf(): Generating inputs from an unmarked
#    frame (e.g. when no newdata) using getDesign
# 3. get_formula: Get formula for submodel type
# 4. get_orig_data(): Get original dataset for use in building model frame
# 5. predict_by_chunk(): Take inputs and generate predictions
# Basic methods are shown below; fit type-specific methods in their own sections

setGeneric("check_predict_arguments", function(object, ...){
  standardGeneric("check_predict_arguments")
})

setMethod("check_predict_arguments", "unmarkedFit",
  function(object, type, newdata, ...){
  # Check if type is supported (i.e., is it in names(object)?)
  check_type(object, type)

  # Check newdata class
  if(!inherits(newdata, c("unmarkedFrame", "data.frame", "RasterLayer", "RasterStack"))){
    stop("newdata must be unmarkedFrame, data.frame, RasterLayer, or RasterStack", call.=FALSE)
  }
  invisible(TRUE)
})

# Check if predict type is valid
check_type <- function(mod, type){
  opts <- names(mod)
  if(type %in% opts) return(invisible(TRUE))
  stop("Valid types are ", paste(opts, collapse=", "), call.=FALSE)
}

# Get X and offset when newdata is umf
setGeneric("predict_inputs_from_umf", function(object, ...){
  standardGeneric("predict_inputs_from_umf")
})

setMethod("predict_inputs_from_umf", "unmarkedFit",
  function(object, type, newdata, na.rm, re.form){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  if(type == "state") list_els <- c("X","Z_state","X.offset")
  if(type == "det") list_els <- c("V","Z_det","V.offset")

  X <- designMats[[list_els[1]]]
  if(is.null(re.form)) X <- cbind(X, designMats[[list_els[2]]])
  offset <- designMats[[list_els[3]]]

  list(X=X, offset=offset)
})

# Get correct individual formula based on type
setGeneric("get_formula", function(object, type, ...){
  standardGeneric("get_formula")
})

setMethod("get_formula", "unmarkedFit", function(object, type, ...){
  if(type == "state"){
    return(as.formula(paste("~", object@formula[3], sep="")))
  } else if(type == "det"){
    return(as.formula(object@formula[[2]]))
  }
  NULL
})

# When newdata is data.frame/raster, get original dataset
# For use in building correct model frame
setGeneric("get_orig_data", function(object, type, ...){
  standardGeneric("get_orig_data")
})

# Note that by default, final year of yearlySiteCov data at each site is dropped
# Because transition probabilities are not estimated for final year
# this is appropriate for dynamic models but not temporary emigration models
# for which the drop_final should be FALSE
setMethod("get_orig_data", "unmarkedFit", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=TRUE)
  datatype <- switch(type, state='site_covs', det='obs_covs')
  clean_covs[[datatype]]
})

# Convert NULL data frames to dummy data frames of proper dimension
# Add site covs to yearlysitecovs, ysc to obs covs, etc.
# Drop final year of ysc if necessary
clean_up_covs <- function(object, drop_final=FALSE){
  M <- numSites(object)
  R <- obsNum(object)
  T <- 1
  J <- R
  is_mult <- methods::.hasSlot(object, "numPrimary")
  if(is_mult){
    T <- object@numPrimary
    J <- R/T
  }

  sc <- siteCovs(object)
  if(is.null(sc)) sc <- data.frame(.dummy=rep(1,M))
  out <- list(site_covs=sc)

  if(is_mult){
    ysc <- yearlySiteCovs(object)
    if(is.null(ysc)) ysc <- data.frame(.dummy2=rep(1,M*T))
    ysc <- cbind(ysc, sc[rep(1:M, each=T),,drop=FALSE])
  }

  if(methods::.hasSlot(object, "obsCovs")){
    oc <- obsCovs(object)
    if(is.null(oc)) oc <- data.frame(.dummy3=rep(1,M*T*J))
    if(is_mult){
      oc <- cbind(oc, ysc[rep(1:(M*T), each=J),,drop=FALSE])
    } else {
      oc <- cbind(oc, sc[rep(1:M, each=J),,drop=FALSE])
    }
    out$obs_covs=oc
  }

  if(is_mult){
    if(drop_final & (T > 1)){
      # Drop final year of data at each site
      # Also drop factor levels only found in last year of data
      ysc <- drop_final_year(ysc, M, T)
    }
    out$yearly_site_covs <- ysc
  }

  out
}

#Remove data in final year of yearlySiteCovs (replacing with NAs)
#then drop factor levels found only in that year
drop_final_year <- function(dat, nsites, nprimary){
  dat[seq(nprimary, nsites*nprimary, by=nprimary), ] <- NA
  dat <- lapply(dat, function(x) x[,drop = TRUE])
  as.data.frame(dat)
}


# Take inputs (most importantly model matrix and offsets) and generate prediction
# done in chunks for speed, 70 was optimal after tests
setGeneric("predict_by_chunk", function(object, ...){
  standardGeneric("predict_by_chunk")
})

setMethod("predict_by_chunk", "unmarkedFit",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){

  if(is.vector(xmat)) xmat <- matrix(xmat, nrow=1)
  nr <- nrow(xmat)
  ind <- rep(1:ceiling(nr/chunk_size), each=chunk_size, length.out=nr)

  # should find a way to keep xmat sparse, but it doesn't
  # work with linearComb
  #x_chunk <- lapply(split(as.data.frame(xmat), ind), as.matrix)
  x_chunk <- lapply(unique(ind),
                    function(i) as.matrix(xmat[ind==i,,drop=FALSE]))

  if(is.null(offsets)) offsets <- rep(0, nr)
  off_chunk <- split(offsets, ind)
  out <- mapply(function(x_i, off_i){
    has_na <- apply(x_i, 1, function(x_i) any(is.na(x_i)))
    # Work around linearComb bug where there can't be NAs in inputs
    x_i[has_na,] <- 0
    off_i[has_na] <- 0
    lc <- linearComb(object, x_i, type, offset=off_i, re.form=re.form)
    if(backTransform) lc <- backTransform(lc)
    out <- data.frame(Predicted=coef(lc))
    if(!is.null(level)){
      se <- SE(lc)
      ci <- confint(lc, level=level)
      out$SE <- se
      out$lower <- ci[,1]
      out$upper <- ci[,2]
    }
    out[has_na,] <- NA
    out
    }, x_chunk, off_chunk, SIMPLIFY=FALSE)
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
})


# Raster handling functions----------------------------------------------------

# Convert a raster into a data frame to use as newdata
newdata_from_raster <- function(rst, vars){
  nd <- raster::as.data.frame(rst)
  # Handle factor rasters
  is_fac <- raster::is.factor(rst)
  rem_string <- paste(paste0("^",names(rst),"_"), collapse="|")
  names(nd)[is_fac] <- gsub(rem_string, "", names(nd)[is_fac])
  # Check if variables are missing
  no_match <- vars[! vars %in% names(nd)]
  if(length(no_match) > 0){
    stop(paste0("Variable(s) ",paste(no_match, collapse=", "), " not found in raster stack"),
         call.=FALSE)
  }
  return(nd)
}

# Convert predict output into a raster
raster_from_predict <- function(pr, orig_rst, appendData){
  new_rast <- data.frame(raster::coordinates(orig_rst), pr)
  new_rast <- raster::stack(raster::rasterFromXYZ(new_rast))
  raster::crs(new_rast) <- raster::crs(orig_rst)
  if(appendData) new_rast <- raster::stack(new_rast, orig_rst)
  new_rast
}


# pcount methods---------------------------------------------------------------

setMethod("check_predict_arguments", "unmarkedFitPCount",
  function(object, type, newdata, ...){
  if(type %in% c("psi", "alpha")){
    stop(paste0(type, " is scalar. Use backTransform instead."), call.=FALSE)
  }
  methods::callNextMethod(object, type, newdata)
})

# Special predict approach for ZIP distribution in pcount
# All other distributions use default method
setMethod("predict_by_chunk", "unmarkedFitPCount",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "state" & object@mixture == "ZIP"){

    out <- data.frame(matrix(NA, nrow(xmat), 4))
    names(out) <- c("Predicted", "SE", "lower", "upper")

    psi.hat <- plogis(coef(object, type="psi"))
    lamEst <- object["state"]
    psiEst <- object["psi"]
    fixedOnly <- !is.null(re.form)
    lam.mle <- coef(lamEst, fixedOnly=fixedOnly)
    lam_vcov <- vcov(lamEst, fixedOnly=fixedOnly)
    if(is.null(offsets)) offsets <- rep(0, nrow(xmat))

    for(i in 1:nrow(xmat)) {
      if(nrow(xmat) > 5000) {
          if(i %% 1000 == 0)
              cat("  doing row", i, "of", nrow(xmat), "\n")
      }
      if(any(is.na(xmat[i,]))) next
      ## for the ZIP model the predicted values on the log scale have us
      ## add log(1-psi.hat) to the normal linear prediction
      out$Predicted[i] <-   xmat[i,] %*% lam.mle + offsets[i] + log(1 - psi.hat)
      ## to compute the approximate SE, I compute the variance of the usual
      ## linear part -- that is easy, and to that I add the variance of
      ## log(1-psi.hat) obtained by the delta approximation
      logit.psi<-coef(object,type="psi")
      #  To do that I took derivative of log(1-psi.hat) using application
      #  of chain rule.... hopefully correctly.
      delta.approx.2ndpart<-   ( ((1/(1-psi.hat))*(exp(logit.psi)/((1+exp(logit.psi))^2)))^2 ) * (SE(psiEst)^2)
      ## now the SE is the sqrt of the whole thing
      out$SE[i]<- sqrt( t(xmat[i,])%*% lam_vcov %*%xmat[i,] + delta.approx.2ndpart   )

      #From Mike Meredith
      alf <- (1 - level) / 2
      crit<-qnorm(c(alf, 1 - alf))
      ci <- out$Predicted[i] + crit * out$SE[i]
      out$lower[i]<- ci[1]
      out$upper[i]<- ci[2]
      if(backTransform){
        out$Predicted[i] <- exp(out$Predicted[i])
        ### If back-transform, delta approx says var = (exp(linear.predictor)^2)*Var(linear.predictor)
        ### also I exponentiate the confidence interval.....
        out$SE[i]<- out$Predicted[i]*out$SE[i]
        ci<-exp(ci)
      }
      out$lower[i] <- ci[1]
      out$upper[i] <- ci[2]
    }
    return(out)
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})


# colext methods---------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitColExt",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  list_el <- switch(type, psi="W", col="X.gam", ext="X.eps", det="V")
  # colext doesn't support offsets
  list(X=designMats[[list_el]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitColExt", function(object, type, ...){
  switch(type, psi=object@psiformula, col=object@gamformula,
                     ext=object@epsformula, det=object@detformula)
})

setMethod("get_orig_data", "unmarkedFitColExt", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=TRUE)
  datatype <- switch(type, psi='site_covs', col='yearly_site_covs',
                     ext='yearly_site_covs', det='obs_covs')
  clean_covs[[datatype]]
})


# occuFP methods---------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitOccuFP",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@detformula, object@FPformula,
                          object@Bformula, object@stateformula, na.rm=na.rm)
  X_idx <- switch(type, state="X", det="V", fp="U", b="W")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=designMats[[off_idx]])
})

setMethod("get_formula", "unmarkedFitOccuFP", function(object, type, ...){
  switch(type, state=object@stateformula, det=object@detformula,
         b=object@Bformula, fp=object@FPformula)
})

setMethod("get_orig_data", "unmarkedFitOccuFP", function(object, type, ...){
  # Get obs data if fp, b, or det
  new_type <- ifelse(type %in% c("fp", "b"), "det", type)
  methods::callNextMethod(object, new_type, ...)
})


# Dail-Madsen model methods----------------------------------------------------

# Includes unmarkedFitPCO, unmarkedFitMMO, unmarkedFitDSO

setMethod("check_predict_arguments", "unmarkedFitDailMadsen",
  function(object, type, newdata, ...){
  if(type %in% c("psi", "alpha", "scale")){
    stop(paste0(type, " is scalar. Use backTransform instead."), call.=FALSE)
  }
  dynamics <- object@dynamics
  immigration <- tryCatch(object@immigration, error=function(e) FALSE)
  if(identical(dynamics, "notrend") & identical(type, "gamma"))
    stop("gamma is a derived parameter for this model: (1-omega)*lambda")
  if(identical(dynamics, "trend") && identical(type, "omega"))
    stop("omega is not a parameter in the dynamics='trend' model")
  if(!immigration && identical(type, "iota"))
    stop("iota is not a parameter in the immigration=FALSE model")
  methods::callNextMethod(object, type, newdata)
})

setMethod("predict_inputs_from_umf", "unmarkedFitDailMadsen",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm=na.rm)
  X_idx <- switch(type, lambda="Xlam", gamma="Xgam", omega="Xom",
                  iota="Xiota", det="Xp")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=designMats[[off_idx]])
})

setMethod("get_formula", "unmarkedFitDailMadsen", function(object, type, ...){
  fl <- object@formlist
  switch(type, lambda=fl$lambdaformula, gamma=fl$gammaformula,
         omega=fl$omegaformula, iota=fl$iotaformula, det=fl$pformula)
})

setMethod("get_orig_data", "unmarkedFitDailMadsen", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=TRUE)
  datatype <- switch(type, lambda='site_covs', gamma='yearly_site_covs',
                     omega='yearly_site_covs', iota='yearly_site_covs',
                     det='obs_covs')
  clean_covs[[datatype]]
})

# This method differs for DSO
setMethod("get_orig_data", "unmarkedFitDSO", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=TRUE)
  datatype <- switch(type, lambda='site_covs', gamma='yearly_site_covs',
                     omega='yearly_site_covs', iota='yearly_site_covs',
                     det='yearly_site_covs')
  clean_covs[[datatype]]
})

# Special handling for ZIP distribution
setMethod("predict_by_chunk", "unmarkedFitDailMadsen",
  function(object, type, level, xmat, offsets, chunk_size, backTransform=TRUE,
           re.form=NULL, ...){
  if(type == "lambda" & object@mixture == "ZIP"){
    warning("Method to compute SE for ZIP model has not been written", call.=FALSE)
    out <- data.frame(matrix(NA, nrow(xmat), 4))
    names(out) <- c("Predicted", "SE", "lower", "upper")
    lam.mle <- coef(object, type="lambda")
    psi.hat <- plogis(coef(object, type="psi"))
    if(is.null(offsets)) offsets <- rep(0, nrow(xmat))
    out$Predicted <- as.numeric(xmat %*% lam.mle + offsets + log(1 - psi.hat))
    if(backTransform) out$Predicted <- exp(out$Predicted)
    return(out)
  }
  methods::callNextMethod(object, type, level, xmat, offsets, chunk_size,
                          backTransform, re.form, ...)
})


# Temporary emigration models--------------------------------------------------

# All inherit from GMM so only one set of methods is required
# (except GDR which has its own predict method right now)

setMethod("predict_inputs_from_umf", "unmarkedFitGMM",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm=na.rm)
  X_idx <- switch(type, lambda="Xlam", phi="Xphi", det="Xdet")
  off_idx <- paste0(X_idx, ".offset")
  list(X=designMats[[X_idx]], offset=designMats[[off_idx]])
})

setMethod("get_formula", "unmarkedFitGMM", function(object, type, ...){
  fl <- object@formlist
  switch(type, lambda=fl$lambdaformula, phi=fl$phiformula, det=fl$pformula)
})

setMethod("get_orig_data", "unmarkedFitGMM", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=FALSE)
  datatype <- switch(type, lambda='site_covs', phi='yearly_site_covs',
                     det='obs_covs')
  clean_covs[[datatype]]
})


# occuTTD----------------------------------------------------------------------

# Identical to colext

setMethod("predict_inputs_from_umf", "unmarkedFitOccuTTD",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  list_el <- switch(type, psi="W", col="X.gam", ext="X.eps", det="V")
  list(X=designMats[[list_el]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitOccuTTD", function(object, type, ...){
  switch(type, psi=object@psiformula, col=object@gamformula,
                     ext=object@epsformula, det=object@detformula)
})

setMethod("get_orig_data", "unmarkedFitOccuTTD", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=TRUE)
  datatype <- switch(type, psi='site_covs', col='yearly_site_covs',
                     ext='yearly_site_covs', det='obs_covs')
  clean_covs[[datatype]]
})


# nmixTTD----------------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitNmixTTD",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formula, na.rm = na.rm)
  list_el <- switch(type, state="W", det="V")
  list(X=designMats[[list_el]], offset=NULL)
})

setMethod("get_formula", "unmarkedFitNmixTTD", function(object, type, ...){
  switch(type, state=object@stateformula, det=object@detformula)
})

setMethod("get_orig_data", "unmarkedFitNmixTTD", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=FALSE)
  datatype <- switch(type, state='site_covs', det='obs_covs')
  clean_covs[[datatype]]
})


# gdistremoval-----------------------------------------------------------------

setMethod("predict_inputs_from_umf", "unmarkedFitGDR",
  function(object, type, newdata, na.rm, re.form=NA){
  designMats <- getDesign(newdata, object@formlist)
  if(type == "lambda") list_els <- c("Xlam","Zlam")
  if(type == "phi") list_els <- c("Xphi","Zphi")
  if(type == "dist") list_els <- c("Xdist","Zdist")
  if(type == "rem") list_els <- c("Xrem", "Zrem")
  X <- designMats[[list_els[1]]]
  if(is.null(re.form)) X <- cbind(X, designMats[[list_els[2]]])
  list(X=X, offset=NULL)
})

setMethod("get_formula", "unmarkedFitGDR", function(object, type, ...){
  fl <- object@formlist
  switch(type, lambda=fl$lambdaformula, phi=fl$phiformula,
         dist=fl$distanceformula, rem=fl$removalformula)
})

setMethod("get_orig_data", "unmarkedFitGDR", function(object, type, ...){
  clean_covs <- clean_up_covs(object@data, drop_final=FALSE)
  datatype <- switch(type, lambda='site_covs', phi='yearly_site_covs',
                     dist='yearly_site_covs', rem='obs_covs')
  clean_covs[[datatype]]
})


# occuMulti--------------------------------------------------------------------

# bespoke predict method since it has numerious unusual options
# and requires bootstrapping

# This method is used by simulate but not by predict
setMethod("get_formula", "unmarkedFitOccuMulti", function(object, type, ...){
  switch(type, state=object@stateformulas,
               det=object@detformulas)
})

setMethod("predict", "unmarkedFitOccuMulti",
     function(object, type, newdata,
              #backTransform = TRUE, na.rm = TRUE,
              #appendData = FALSE,
              se.fit=TRUE, level=0.95, species=NULL, cond=NULL, nsims=100,
              ...)
  {

  type <- match.arg(type, c("state", "det"))

  if(is.null(hessian(object))){
    se.fit = FALSE
  }

  species <- name_to_ind(species, names(object@data@ylist))
  cond <- name_to_ind(cond, names(object@data@ylist))

  if(missing(newdata)){
    newdata <- NULL
  } else {
    if(! class(newdata) %in% c('data.frame')){
      stop("newdata must be a data frame")
    }
  }

  maxOrder <- object@call$maxOrder
  if(is.null(maxOrder)) maxOrder <- length(object@data@ylist)
  dm <- getDesign(object@data,object@detformulas,object@stateformulas,
                  maxOrder, na.rm=F, newdata=newdata, type=type)

  params <- coef(object)
  low_bound <- (1-level)/2
  up_bound <- level + (1-level)/2


  if(type=="state"){
    N <- nrow(dm$dmOcc[[1]]); nF <- dm$nF; dmOcc <- dm$dmOcc;
    fStart <- dm$fStart; fStop <- dm$fStop; fixed0 <- dm$fixed0
    t_dmF <- t(dm$dmF)

    calc_psi <- function(params){

      f <- matrix(NA,nrow=N,ncol=nF)
      index <- 1
      for (i in 1:nF){
        if(fixed0[i]){
          f[,i] <- 0
        } else {
          f[,i] <- dmOcc[[index]] %*% params[fStart[index]:fStop[index]]
          index <- index + 1
        }
      }
      psi <- exp(f %*% t_dmF)
      as.matrix(psi/rowSums(psi))
    }

    psi_est <- calc_psi(params)

    if(se.fit){
      message('Bootstrapping confidence intervals with ',nsims,' samples')
      Sigma <- vcov(object)
      samp <- array(NA,c(dim(psi_est),nsims))
      for (i in 1:nsims){
        samp[,,i] <- calc_psi(mvrnorm(1, coef(object), Sigma))
      }
    }

    if(!is.null(species)){

      sel_col <- species

      if(!is.null(cond)){
        if(any(sel_col %in% abs(cond))){
          stop("Species can't be conditional on itself")
        }
        ftemp <- object@data@fDesign
        swap <- -1*cond[which(cond<0)]
        ftemp[,swap] <- 1 - ftemp[,swap]
        num_inds <- apply(ftemp[,c(sel_col,abs(cond))] == 1,1,all)
        denom_inds <- apply(ftemp[,abs(cond),drop=F] == 1,1,all)
        est <- rowSums(psi_est[,num_inds,drop=F]) /
          rowSums(psi_est[,denom_inds, drop=F])
        if(se.fit){
          samp_num <- apply(samp[,num_inds,,drop=F],3,rowSums)
          samp_denom <- apply(samp[,denom_inds,,drop=F],3,rowSums)
          samp <- samp_num / samp_denom
        }

      } else {
        num_inds <- apply(object@data@fDesign[,sel_col,drop=FALSE] == 1,1,all)
        est <- rowSums(psi_est[,num_inds,drop=F])
        if(se.fit){
          samp <- samp[,num_inds,,drop=F]
          samp <- apply(samp, 3, rowSums)
        }
      }

      if(se.fit){
        if(!is.matrix(samp)) samp <- matrix(samp, nrow=1)
        boot_se <- apply(samp,1,sd, na.rm=T)
        boot_low <- apply(samp,1,quantile,low_bound, na.rm=T)
        boot_up <- apply(samp,1,quantile,up_bound, na.rm=T)
      } else{
        boot_se <- boot_low <- boot_up <- NA
      }
      return(data.frame(Predicted=est,
                        SE=boot_se,
                        lower=boot_low,
                        upper=boot_up))

    } else {
      codes <- apply(dm$z,1,function(x) paste(x,collapse=""))
      colnames(psi_est)  <- paste('psi[',codes,']',sep='')
      if(se.fit){
        boot_se <- apply(samp,c(1,2),sd, na.rm=T)
        boot_low <- apply(samp,c(1,2),quantile,low_bound, na.rm=T)
        boot_up <- apply(samp,c(1,2),quantile,up_bound, na.rm=T)
        colnames(boot_se) <- colnames(boot_low) <- colnames(boot_up) <-
          colnames(psi_est)
      } else {
        boot_se <- boot_low <- boot_up <- NA
      }
      return(list(Predicted=psi_est,
                  SE=boot_se,
                  lower=boot_low,
                  upper=boot_up))
    }
  }

  if(type=="det"){
    S <- dm$S; dmDet <- dm$dmDet
    dStart <- dm$dStart; dStop <- dm$dStop

    out <- list()
    for (i in 1:S){
      #Subset estimate to species i
      inds <- dStart[i]:dStop[i]
      new_est <- object@estimates@estimates$det
      new_est@estimates <- coef(object)[inds]
      new_est@fixed <- 1:length(inds)
      if(se.fit){
        new_est@covMat <- vcov(object)[inds,inds,drop=FALSE]
        new_est@covMatBS <- object@covMatBS[inds,inds,drop=FALSE]
      } else{
        new_est@covMat <- matrix(NA, nrow=length(inds), ncol=length(inds))
        new_est@covMatBS <- matrix(NA, nrow=length(inds), ncol=length(inds))
      }

      prmat <- t(apply(dmDet[[i]], 1, function(x){
                    bt <- backTransform(linearComb(new_est, x))
                    if(!se.fit){
                      return(c(Predicted=bt@estimate, SE=NA, lower=NA, upper=NA))
                    }
                    ci <- confint(bt, level=level)
                    names(ci) <- c("lower", "upper")
                    c(Predicted=bt@estimate, SE=SE(bt), ci)
                  }))
      rownames(prmat) <- NULL
      out[[i]] <- as.data.frame(prmat)
    }
    names(out) <- names(object@data@ylist)
    if(!is.null(species)){
      return(out[[species]])
    }
    return(out)
  }
  stop("type must be 'det' or 'state'")
})


# occuMS-----------------------------------------------------------------------

# bespoke predict method since it requires bootstrapping

# This method is used by simulate by not by predict
setMethod("get_formula", "unmarkedFitOccuMS", function(object, type, ...){
  switch(type, psi=object@psiformulas, phi=object@phiformulas,
               det=object@detformulas)
})

setMethod("predict", "unmarkedFitOccuMS",
     function(object, type, newdata,
              #backTransform = TRUE, na.rm = TRUE,
              #appendData = FALSE,
              se.fit=TRUE, level=0.95, nsims=100, ...)
{

  #Process input---------------------------------------------------------------
  if(! type %in% c("psi","phi", "det")){
    stop("type must be 'psi', 'phi', or 'det'")
  }

  if(is.null(hessian(object))){
    se.fit = FALSE
  }

  if(missing(newdata)){
    newdata <- NULL
  } else {
    if(! class(newdata) %in% c('data.frame')){
      stop("newdata must be a data frame")
    }
  }

  S <- object@data@numStates
  gd <- getDesign(object@data,object@psiformulas,object@phiformulas,
                  object@detformulas, object@parameterization, na.rm=F,
                  newdata=newdata, type=type)

  #Index guide used to organize p values
  guide <- matrix(NA,nrow=S,ncol=S)
  guide <- lower.tri(guide,diag=TRUE)
  guide[,1] <- FALSE
  guide <- which(guide,arr.ind=TRUE)
  #----------------------------------------------------------------------------

  #Utility functions-----------------------------------------------------------
  #Get matrix of linear predictor values
  get_lp <- function(params, dm_list, ind){
    L <- length(dm_list)
    out <- matrix(NA,nrow(dm_list[[1]]),L)
    for (i in 1:L){
      out[,i] <- dm_list[[i]] %*% params[ind[i,1]:ind[i,2]]
    }
    out
  }

  #Get SE/CIs for conditional binomial using delta method
  split_estimate <- function(object, estimate, inds, se.fit){
    out <- estimate
    out@estimates <- coef(object)[inds]
    if(se.fit){
      out@covMat <- vcov(object)[inds,inds,drop=FALSE]
    } else{
      out@covMat <- matrix(NA, nrow=length(inds), ncol=length(inds))
    }
    out
  }

  lc_to_predict <- function(object, estimate, inds, dm, level, se.fit){

    new_est <- split_estimate(object, estimate, inds[1]:inds[2], se.fit)

    out <- t(apply(dm, 1, function(x){
      bt <- backTransform(linearComb(new_est, x))
      if(!se.fit) return(c(Predicted=bt@estimate, SE=NA, lower=NA, upper=NA))
      ci <- confint(bt, level=level)
      names(ci) <- c("lower", "upper")
      c(Predicted=bt@estimate, SE=SE(bt), ci)
    }))
    rownames(out) <- NULL
    as.data.frame(out)
  }


  #Calculate row-wise multinomial logit prob
  #implemented in C++ below as it is quite slow
  get_mlogit_R <- function(lp_mat){
    if(type == 'psi'){
      out <- cbind(1,exp(lp_mat))
      out <- out/rowSums(out)
      out <- out[,-1]
    } else if(type == 'phi'){ #doesn't work
      np <- nrow(lp_mat)
      out <- matrix(NA,np,ncol(lp_mat))
      ins <- outer(1:S, 1:S, function(i,j) i!=j)
      for (i in 1:np){
        phimat <- diag(S)
        phimat[ins] <- exp(lp_mat[i,])
        phimat <- t(phimat)
        phimat <- phimat/rowSums(phimat)
        out[i,] <- phimat[ins]
      }
    } else {
      R <- nrow(lp_mat)
      out <- matrix(NA,R,ncol(lp_mat))
      for (i in 1:R){
        sdp <- matrix(0,nrow=S,ncol=S)
        sdp[guide] <- exp(lp_mat[i,])
        sdp[,1] <- 1
        sdp <- sdp/rowSums(sdp)
        out[i,] <- sdp[guide]
      }
    }
    out
  }

  get_mlogit <- function(lp_mat){
    .Call("get_mlogit",
         lp_mat, type, S, guide-1)
  }

  #----------------------------------------------------------------------------

  if(type=="psi"){
    dm_list <- gd$dm_state
    ind <- gd$state_ind
    est <- object@estimates@estimates$state
  } else if(type=="phi"){
    dm_list <- gd$dm_phi
    ind <- gd$phi_ind
    est <- object@estimates@estimates$transition
  } else {
    dm_list <- gd$dm_det
    ind <- gd$det_ind
    est <- object@estimates@estimates$det
  }

  P <- length(dm_list)

  low_bound <- (1-level)/2
  z <- qnorm(low_bound,lower.tail=F)

  out <- vector("list", P)
  names(out) <- names(dm_list)

  if(object@parameterization == 'condbinom'){
    out <- lapply(1:length(dm_list), function(i){
      lc_to_predict(object, est, ind[i,], dm_list[[i]], level, se.fit)
    })
    names(out) <- names(dm_list)
    return(out)

  } else if (object@parameterization == "multinomial"){
    lp <- get_lp(coef(object), dm_list, ind)
    pred <- get_mlogit(lp)

    M <- nrow(pred)
    upr <- lwr <- se <- matrix(NA,M,P)

    if(se.fit){
      message('Bootstrapping confidence intervals with',nsims,'samples')

      sig <- vcov(object)
      param_mean <- coef(object)
      rparam <- mvrnorm(nsims, param_mean, sig)

      get_pr <- function(i){
        lp <- get_lp(rparam[i,], dm_list, ind)
        get_mlogit(lp)
      }
      samp <- sapply(1:nsims, get_pr, simplify='array')

      for (i in 1:M){
        for (j in 1:P){
          dat <- samp[i,j,]
          se[i,j] <- sd(dat, na.rm=TRUE)
          quants <- quantile(dat, c(low_bound, (1-low_bound)),na.rm=TRUE)
          lwr[i,j] <- quants[1]
          upr[i,j] <- quants[2]
        }
      }

    }
  }

  for (i in 1:P){
    out[[i]] <- data.frame(Predicted=pred[,i], SE=se[,i],
                           lower=lwr[,i], upper=upr[,i])
  }

  out
})
