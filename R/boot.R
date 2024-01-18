

# ----------------------- PARAMETRIC BOOTSTRAP --------------------------

setGeneric("parboot",
           def = function(object, ...) {
             standardGeneric("parboot")
           })


setClass("parboot",
         representation(call = "call",
                        t0 = "numeric",
                        t.star = "matrix"))

setGeneric("replaceY", function(object, newY, replNA = TRUE, ...){
           standardGeneric("replaceY")})
setMethod("replaceY", "unmarkedFrame", function(object, newY, replNA=TRUE, ...){
    if(replNA) is.na(newY) <- is.na(object@y)
    object@y <- newY
    object
})
setMethod("replaceY", "unmarkedFrameOccuMulti",
          function(object, newY, replNA=TRUE, ...){
      if(replNA){
        newY <- mapply(function(x, y){ is.na(x) <- is.na(y); x},
                       newY , object@ylist, SIMPLIFY=FALSE)
      }
      object@ylist <- newY
      object
})


setMethod("parboot", "unmarkedFit",  function(object, statistic=SSE, nsim=10,
          report, seed = NULL, parallel = FALSE, ncores, ...){

  if(!missing(report)){
    warning("report argument is non-functional and will be deprecated in the next version", call.=FALSE)
  }

  dots <- list(...)
  call <- match.call(call = sys.call(-1))
  stopifnot(is.function(statistic))
  starts <- as.numeric(coef(object))
  # Get rid of starting values if model was fit with TMB
  if(methods::.hasSlot(object, "TMB") && !is.null(object@TMB)) starts <- NULL

  t0 <- statistic(object, ...)

  simList <- simulate(object, nsim = nsim, na.rm = FALSE)

  availcores <- parallel::detectCores() - 1
  if(missing(ncores) || ncores > availcores) ncores <- availcores

  cl <- NULL
  if(parallel){
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, library(unmarked))
    env_vars <- c("dots", "replaceY")
    fm.nms <- all.names(object@call)
    if (!any(grepl("~", fm.nms))) env_vars <- c(env_vars, fm.nms[2])
    if(.hasSlot(object@data, "piFun")) env_vars <- c(env_vars, object@data@piFun)
    parallel::clusterExport(cl, env_vars, envir = environment())
    parallel::clusterEvalQ(cl, list2env(dots))
  }

  run_sim <- function(x, object, statistic, starts, t0, ...){
    simdata <- replaceY(object@data, x)
    tryCatch({
      #if(runif(1,0,1) < 0.5) stop("fail") # for testing error trapping
      fit <- update(object, data = simdata, starts = starts, se = FALSE)
      statistic(fit, ...)
    }, error=function(e){
      t0[] <- NA
      t0
    })
  }

  # Uses pbapply if available, or parSapply if not (see utils.R)
  t.star <- t(sapply2(simList, run_sim, object=object,
                              statistic=statistic, starts=starts, t0=t0,
                              cl=cl, ...))
  if(length(t0) == 1) t.star <- matrix(t.star, ncol=1)

  failed <- apply(t.star, 1, function(x) any(is.na(x)))
  if(sum(failed) > 0){
    warning(paste0("Model fitting failed in ",sum(failed), " sims."), call.=FALSE)
    t.star <- t.star[!failed,,drop=FALSE]
  }

  new("parboot", call = call, t0 = t0, t.star = t.star)

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
    cat("\nCall:", deparse(object@call, width.cutoff=500), fill=T)
    cat("\nParametric Bootstrap Statistics:\n")
    print(stats, digits=3)
    cat("\nt_B quantiles:\n")
    print(t(apply(t.star, 2, quantile,
        probs=c(0, 2.5, 25, 50, 75, 97.5, 100) / 100)), digits=2)
    cat("\nt0 = Original statistic computed from data\n")
    cat("t_B = Vector of bootstrap samples\n\n")
})




setMethod("plot", signature(x="parboot", y="missing"),
    function(x, y, xlab, main = "Parametric Bootstrapped Samples", xlim,
        ...)
{
    t.star <- x@t.star
    t0 <- x@t0
    if(length(t0) > 1) {
      oldask <- devAskNewPage(ask = dev.interactive(orNone = TRUE))
      on.exit(devAskNewPage(oldask))
    }
    for(i in 1:length(t0)) {
        if(missing(xlab))
            xlab <- colnames(t.star)[i]
        h <- hist(t.star[,i], plot = FALSE)
        if(missing(xlim))
            xl <- c(min(h$breaks[1], t0[i]), max(max(h$breaks), t0[i]))
        else
            xl <- xlim
        hist(t.star[,i], xlab=xlab, xlim = xl, main = main, ...)
        abline(v=t0[i], lty=2)
    }
})


# ----------------------- Nonparametric bootstrapping -------------------

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
    designMats <- getDesign(data, formula) # bootstrap after removing sites
    removed.sites <- designMats$removed.sites
    if(length(removed.sites)>0)
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
            obs.per.site <- lapply(1:nrow(y), function(i) which(!is.na(y[i,])))
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
                                 replicate(B, boot.iter(),
                                           simplify = FALSE))
    coefs <- t(sapply(object@bootstrapSamples,
                      function(x) coef(x)))
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
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="both")
})


setMethod("nonparboot", "unmarkedFitPCount",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})


setMethod("nonparboot", "unmarkedFitMPois",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})


setMethod("nonparboot", "unmarkedFitDS",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})


setMethod("nonparboot", "unmarkedFitOccuRN",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="both")
})



setMethod("nonparboot", "unmarkedFitGMM",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})

setMethod("nonparboot", "unmarkedFitGDS",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})




setMethod("nonparboot", "unmarkedFitDailMadsen",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})





setMethod("nonparboot", "unmarkedFitColExt",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
#    browser()
    if (identical(B, 0) && !is.null(object@bootstrapSamples))
        return(object)
    if (B <= 0 && is.null(object@bootstrapSamples))
        stop("B must be greater than 0 when fit has no bootstrap samples.")
    data <- object@data
    psiParms <- coef(object, 'psi')
    detParms <- coef(object, 'det')
    colParms <- coef(object, 'col')
    extParms <- coef(object, 'ext')

    # bootstrap only after removing sites
    designMats <- getDesign(object@data, formula=object@formula)
    removed.sites <- designMats$removed.sites
    if(length(removed.sites) > 0) {
        sites <- 1:nrow(getY(data))
        keep <- which(!sites %in% removed.sites)
        data <- data[keep,]
        }
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
    if(!keepOldSamples)
        object@bootstrapSamples <- NULL
    object@bootstrapSamples <- c(object@bootstrapSamples,
        replicate(B, boot.iter(), simplify = FALSE))
    coefs <- t(sapply(object@bootstrapSamples, function(x) coef(x)))
    v <- cov(coefs)
    object@covMatBS <- v
    inds <- .estimateInds(object)
    for(est in names(inds)) {
         v.est <- v[inds[[est]], inds[[est]], drop = FALSE]
         object@estimates@estimates[[est]]@covMatBS <- v.est
         }
    smoothed.occ <- t(sapply(object@bootstrapSamples,
         function(x) x@smoothed.mean[1,]))
    smoothed.unocc <- t(sapply(object@bootstrapSamples,
         function(x) x@smoothed.mean[2,]))
    object@smoothed.mean.bsse <-
         rbind(sqrt(diag(cov(smoothed.occ))),
               sqrt(diag(cov(smoothed.unocc))))
    projected.occ <- t(sapply(object@bootstrapSamples,
         function(x) x@projected.mean[1,]))
    projected.unocc <- t(sapply(object@bootstrapSamples,
         function(x) x@projected.mean[2,]))
             object@projected.mean.bsse <-
                rbind(sqrt(diag(cov(projected.occ))),
                      sqrt(diag(cov(projected.unocc))))
    return(object)
})


setMethod("nonparboot", "unmarkedFitOccuPEN",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
#    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
#                   bsType="site")
    bsType <- "site"
    if (identical(B, 0) && !is.null(object@bootstrapSamples)) {
        return(object)
    }
    if (B <= 0 && is.null(object@bootstrapSamples)) {
        stop("B must be greater than 0 when fit has no bootstrap samples.")
    }
    data <- object@data
    formula <- object@formula
    designMats <- getDesign(data, formula) # bootstrap after removing sites
    removed.sites <- designMats$removed.sites
    if(length(removed.sites)>0)
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
            obs.per.site <- lapply(1:nrow(y), function(i) which(!is.na(y[i,])))
            obs <- lapply(obs.per.site,
                          function(obs) sample(obs, replace = TRUE))
            data.b <- data.b[obs]
        }
        fm <- update(object, data = data.b)
        return(fm)
    }
    if (!keepOldSamples) {
        object@bootstrapSamples <- NULL
    }
    object@bootstrapSamples <- c(object@bootstrapSamples,
                                 replicate(B, boot.iter(),
                                           simplify = FALSE))
    coefs <- t(sapply(object@bootstrapSamples,
                      function(x) coef(x)))
    v <- cov(coefs)
    object@covMatBS <- v
    inds <- .estimateInds(object)
    for (est in names(inds)) {
        v.est <- v[inds[[est]], inds[[est]], drop = FALSE]
        object@estimates@estimates[[est]]@covMatBS <- v.est
    }
    object


})


setMethod("nonparboot", "unmarkedFitOccuPEN_CV",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
#    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
#                   bsType="site")
    bsType <- "site"
    if (identical(B, 0) && !is.null(object@bootstrapSamples)) {
        return(object)
    }
    if (B <= 0 && is.null(object@bootstrapSamples)) {
        stop("B must be greater than 0 when fit has no bootstrap samples.")
    }
    data <- object@data
    formula <- object@formula
    designMats <- getDesign(data, formula) # bootstrap after removing sites
    removed.sites <- designMats$removed.sites
    if(length(removed.sites)>0)
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
            obs.per.site <- lapply(1:nrow(y), function(i) which(!is.na(y[i,])))
            obs <- lapply(obs.per.site,
                          function(obs) sample(obs, replace = TRUE))
            data.b <- data.b[obs]
        }
	if (object@pen.type=="MPLE") {
	  MPLElambda = computeMPLElambda(formula,data.b)
	  fm <- update(object, data = data.b,lambda=MPLElambda)
	} else {
          fm <- update(object, data = data.b)
	}
        return(fm)
    }
    if (!keepOldSamples) {
        object@bootstrapSamples <- NULL
    }
    object@bootstrapSamples <- c(object@bootstrapSamples,
                                 replicate(B, boot.iter(),
                                           simplify = FALSE))
    coefs <- t(sapply(object@bootstrapSamples,
                      function(x) coef(x)))
    v <- cov(coefs)
    object@covMatBS <- v
    inds <- .estimateInds(object)
    for (est in names(inds)) {
        v.est <- v[inds[[est]], inds[[est]], drop = FALSE]
        object@estimates@estimates[[est]]@covMatBS <- v.est
    }
    object


})


setMethod("nonparboot", "unmarkedFitOccuTTD",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})



setMethod("nonparboot", "unmarkedFitOccuMulti",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    bsType <- "site"
    if (identical(B, 0) && !is.null(object@bootstrapSamples)) {
        return(object)
    }
    if (B <= 0 && is.null(object@bootstrapSamples)) {
        stop("B must be greater than 0 when fit has no bootstrap samples.")
    }
    data <- object@data
    M <- numSites(data)
    boot.iter <- function() {
      finish <- FALSE
      while(!finish){
        sites <- sort(sample(1:M, M, replace = TRUE))
        data.b <- data[sites,]
        ran <- TRUE
        tryCatch(fm <- update(object, data = data.b, se=FALSE), error=function(e) ran <<-FALSE)
        if(!ran) next
        finish <- fm@opt$convergence == 0
      }
      return(fm)
    }
    if (!keepOldSamples) {
        object@bootstrapSamples <- NULL
    }
    object@bootstrapSamples <- c(object@bootstrapSamples,
                                 replicate(B, boot.iter(),
                                           simplify = FALSE))
    coefs <- t(sapply(object@bootstrapSamples,
                      function(x) coef(x)))
    v <- cov(coefs)
    object@covMatBS <- v
    inds <- .estimateInds(object)
    for (est in names(inds)) {
        v.est <- v[inds[[est]], inds[[est]], drop = FALSE]
        object@estimates@estimates[[est]]@covMatBS <- v.est
    }
    object
})

setMethod("nonparboot", "unmarkedFitNmixTTD",
    function(object, B = 0, keepOldSamples = TRUE, ...)
{
    callNextMethod(object, B=B, keepOldSamples=keepOldSamples,
                   bsType="site")
})












# ----------------------- Helper functions -------------------------------

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
      return(c(prev.list, list(seq(prev.max+1, prev.max +
                                   estimateLengths[type]))))
    }
  }
  retlist <- estimateInds(length(estimateLengths))
  names(retlist) <- names(umf)
  retlist
}


