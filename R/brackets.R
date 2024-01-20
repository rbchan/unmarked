# i is the vector of sites to extract

# Select a subset of sites i---------------------------------------------------
# i is a numeric vector of sites

setMethod("[", c("unmarkedFrame", "numeric", "missing", "missing"),
  function(x, i, j){
  select_sites_internal(x, i)
})

setGeneric("select_sites_internal", function(x, i){
  standardGeneric("select_sites_internal")
})

setMethod("select_sites_internal", "unmarkedFrame",
  function(x, i){

  M <- numSites(x)
  if(length(i) == 0) return(x)
  if(any(i < 0) && any(i > 0))
      stop("i must be all positive or all negative indices.")
  if(all(i < 0)) { # if i is negative, then convert to positive
      i <- (1:M)[i]
      }
  y <- getY(x)[i,]
  if (length(i) == 1) {
      y <- t(y)
      }
  siteCovs <- siteCovs(x)
  obsCovs <- obsCovs(x)
  if (!is.null(siteCovs)) {
      siteCovs <- siteCovs(x)[i, , drop = FALSE]
      }
  if (!is.null(obsCovs)) {
      R <- obsNum(x)
      .site <- rep(1:M, each = R)
      oc <- lapply(i, function(ind){
       obsCovs[.site==ind,,drop=FALSE]
      })
      obsCovs <- do.call(rbind, oc)
      }
  umf <- x
  umf@y <- y
  umf@siteCovs <- siteCovs
  umf@obsCovs <- obsCovs
  umf
})


setMethod("select_sites_internal", "unmarkedFrameOccuMulti",
    function(x, i)
{
    if(length(i) == 0) return(x)
    M <- numSites(x)

    ylist <- lapply(x@ylist,function(x) x[i,,drop=F])
    siteCovs <- siteCovs(x)
    obsCovs <- obsCovs(x)
    if (!is.null(siteCovs)) {
        siteCovs <- siteCovs(x)[i, , drop = FALSE]
        }
    if (!is.null(obsCovs)) {
        R <- obsNum(x)
        .site <- rep(1:M, each = R)
        oc <- lapply(i, function(ind){
         obsCovs[.site==ind,,drop=FALSE]
        })
        obsCovs <- do.call(rbind, oc)
        }
    umf <- x
    umf@y <- ylist[[1]]
    umf@ylist <- ylist
    umf@siteCovs <- siteCovs
    umf@obsCovs <- obsCovs
    umf
})


setMethod("select_sites_internal", "unmarkedMultFrame",
		function(x, i)
{
    M <- numSites(x)
    if(length(i) == 0) return(x)
    if(any(i < 0) && any(i > 0))
        stop("i must be all positive or all negative indices.")
    if(all(i < 0)) { # if i is negative, then convert to positive
        i <- (1:M)[i]
        }
    oldy <- getY(x)
    y <- oldy[i,]
    siteCovs <- siteCovs(x)
    obsCovs <- obsCovs(x)
    if (!is.null(siteCovs)) {
        siteCovs <- siteCovs(x)[i, , drop = FALSE]
        }
    if (!is.null(obsCovs)) {
        R <- obsNum(x)
        .site <- rep(1:M, each = obsNum(x)) #NULL     ## testing
        oc <- lapply(i, function(ind){
         obsCovs[.site==ind,,drop=FALSE]
        })
        obsCovs <- do.call(rbind, oc)
        }
    u <- unmarkedMultFrame(y=matrix(y, ncol=ncol(oldy)),
                           siteCovs=siteCovs,
                           obsCovs=obsCovs,
                           numPrimary=x@numPrimary)
    ysc <- x@yearlySiteCovs
    if(!is.null(ysc)) {
        T <- x@numPrimary
        sites <- rep(1:M, each=T)
        keep <- as.vector(sapply(i, function(x) which(sites %in% x)))
        ysc <- ysc[keep,, drop=FALSE]
        u@yearlySiteCovs <- ysc
        }
    u

})


setMethod("select_sites_internal", "unmarkedFrameOccuMS",
    function(x, i)
{
  multf <- callNextMethod(x, i)
  unmarkedFrameOccuMS(y=getY(multf), siteCovs=siteCovs(multf),
                      yearlySiteCovs=yearlySiteCovs(multf),
                      obsCovs=obsCovs(multf),
                      numPrimary=x@numPrimary)
})


setMethod("select_sites_internal", "unmarkedFrameGMM",
		function(x, i)
{
    M <- nrow(x@y)
    y <- x@y[i,,drop=FALSE]
    R <- obsNum(x)
    T <- x@numPrimary

    sc <- siteCovs(x)[i,,drop=FALSE]

    ysc_ind <- rep(1:M, each=T)
    ysc <- do.call("rbind", lapply(i, function(ind){
      yearlySiteCovs(x)[ysc_ind == ind,,drop=FALSE]
    }))

    oc_ind <- rep(1:M, each=R)
    oc <- do.call("rbind", lapply(i, function(ind){
      obsCovs(x)[oc_ind == ind,,drop=FALSE]
    }))

    unmarkedFrameGMM(y=y, siteCovs=sc,
                     yearlySiteCovs=ysc,
                     obsCovs=oc,
                     piFun=x@piFun, type=x@samplingMethod,
                     obsToY=x@obsToY, numPrimary=x@numPrimary)
})


setMethod("select_sites_internal", "unmarkedFrameGPC",
		function(x, i)
{
    multf <- callNextMethod(x, i) # unmarkedMultFrame
    class(multf) <- "unmarkedFrameGPC"
    multf
})


setMethod("select_sites_internal", "unmarkedFrameGPC",
		function(x, i)
{
    multf <- as(x, "unmarkedMultFrame")
    out <- callNextMethod(multf, i) # unmarkedMultFrame
    as(out, "unmarkedFrameGPC")
})


setMethod("select_sites_internal", "unmarkedFrameGDS",
		function(x, i)
{
    multf <- callNextMethod(x, i) # unmarkedMultFrame
    sur <- x@survey
    if(sur=="line")
        unmarkedFrameGDS(y=getY(multf), siteCovs=siteCovs(multf),
                         yearlySiteCovs=yearlySiteCovs(multf),
                         numPrimary=x@numPrimary,
                         dist.breaks=x@dist.breaks,
                         tlength=x@tlength[i],
                         survey=sur,
                         unitsIn=x@unitsIn)
    else if(sur=="point")
        unmarkedFrameGDS(y=getY(multf), siteCovs=siteCovs(multf),
                         yearlySiteCovs=yearlySiteCovs(multf),
                         numPrimary=x@numPrimary,
                         dist.breaks=x@dist.breaks,
                         survey=sur,
                         unitsIn=x@unitsIn)
})


setMethod("select_sites_internal", "unmarkedFramePCO",
		function(x, i)
{
    multf <- callNextMethod(x, i) # unmarkedMultFrame
    unmarkedFramePCO(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=x@numPrimary,
                     primaryPeriod=x@primaryPeriod[i,,drop=FALSE])
})


setMethod("select_sites_internal", "unmarkedFrameOccuTTD",
		function(x, i)
{
    multf <- callNextMethod(x, i) # unmarkedMultFrame
    unmarkedFrameOccuTTD(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=x@numPrimary,
                     surveyLength=x@surveyLength[i,,drop=FALSE])
})


setMethod("select_sites_internal", "unmarkedFrameDSO",
		function(x, i)
{
    multf <- callNextMethod(x, i) # unmarkedMultFrame
    sur <- x@survey
    pp <- x@primaryPeriod[i,,drop=FALSE]
    if(sur=="line")
        unmarkedFrameDSO(y=getY(multf), siteCovs=siteCovs(multf),
                         yearlySiteCovs=yearlySiteCovs(multf),
                         numPrimary=x@numPrimary,
                         dist.breaks=x@dist.breaks,
                         tlength=x@tlength[i],
                         survey=sur,
                         unitsIn=x@unitsIn,
                         primaryPeriod=pp)
    else if(sur=="point")
        unmarkedFrameDSO(y=getY(multf), siteCovs=siteCovs(multf),
                         yearlySiteCovs=yearlySiteCovs(multf),
                         numPrimary=x@numPrimary,
                         dist.breaks=x@dist.breaks,
                         survey=sur,
                         unitsIn=x@unitsIn,
                         primaryPeriod=pp)
})





# Select a subset of observations j -------------------------------------------
# j is a numeric vector

### RBC: Why??? this doesn't allow umf[,c(1,1)]
setMethod("[", c("unmarkedFrame", "missing", "numeric", "missing"),
		function(x, i, j){
  select_obs_internal(x, j)
})

setGeneric("select_obs_internal", function(x, j){
  standardGeneric("select_obs_internal")
})


setMethod("select_obs_internal", "unmarkedFrame",
		function(x, j)
{
    y <- getY(x)
    obsCovs <- obsCovs(x)
    obsToY <- obsToY(x)
    obs.remove <- rep(TRUE, obsNum(x))
    obs.remove[j] <- FALSE
    y.remove <- t(obs.remove) %*% obsToY > 0
    y <- y[,!y.remove, drop=FALSE]
    obsCovs <- obsCovs[!rep(obs.remove, numSites(x)),, drop=FALSE]
    x@obsCovs <- obsCovs
    x@y <- y
    x@obsToY <- obsToY[!obs.remove,!y.remove, drop=FALSE]
    x
})


setMethod("select_obs_internal", "unmarkedFrameOccuMulti",
		function(x, j)
{
    y <- getY(x)
    obsCovs <- obsCovs(x)
    obsToY <- obsToY(x)
    obs.remove <- rep(TRUE, obsNum(x))
    obs.remove[j] <- FALSE
    y.remove <- t(obs.remove) %*% obsToY > 0
    ylist <- lapply(x@ylist, function(z) z[,!y.remove, drop=F])
    obsCovs <- obsCovs[!rep(obs.remove, numSites(x)),, drop=FALSE]

    x@obsCovs <- obsCovs
    x@y <- ylist[[1]]
    x@ylist <- ylist
    x@obsToY <- obsToY[!obs.remove,!y.remove, drop=FALSE]
    x
})


## for multframes, must remove years at a time
setMethod("select_obs_internal", "unmarkedMultFrame",
		function(x, j)
{
    J <- obsNum(x)/x@numPrimary
    obs <- rep(1:x@numPrimary, each = J)
    years <- 1:x@numPrimary
    numPrimary <- length(j)
    obsj <- match(obs, j)
    j2 <- which(!is.na(obsj))
    u <- callNextMethod(x, j2)
    ysc <- yearlySiteCovs(x)
    if(!is.null(ysc)) {
        ysc <- ysc[rep(!is.na(match(years, j)), nrow(getY(x))),, drop=FALSE]
        u@yearlySiteCovs <- ysc
        }
    u@numPrimary <- numPrimary
    return(u)
})


setMethod("select_obs_internal", "unmarkedFramePCO",
		function(x, j)
{
    multf <- callNextMethod(x, j) # unmarkedMultFrame
    unmarkedFramePCO(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=length(j),
                     primaryPeriod=x@primaryPeriod[,j,drop=FALSE])
})




setMethod("select_obs_internal", "unmarkedFrameOccuTTD",
		function(x, j)
{

    if(any(j>x@numPrimary)) stop("Can't select primary periods that don't exist", call.=FALSE)
    if(!all(j>0)) stop("All indices must be positive", call.=FALSE)
    
    R <- ncol(getY(x))/x@numPrimary
    pp_vec <- rep(1:x@numPrimary, each=R)
    keep_cols <- which(pp_vec%in%j)
    y <- getY(x)[,keep_cols,drop=FALSE]
    sl <- x@surveyLength[,keep_cols,drop=FALSE]

    pp_vec2 <- rep(1:x@numPrimary, numSites(x))
    keep_rows <- which(pp_vec2 %in% j)
    ysc <- yearlySiteCovs(x)[keep_rows,,drop=FALSE]

    obs_vec <- rep(rep(1:x@numPrimary, each = R), numSites(x))
    keep_rows <- which(obs_vec %in% j)
    oc <- obsCovs(x)[keep_rows,,drop=FALSE]

    unmarkedFrameOccuTTD(y=y, surveyLength=sl, siteCovs=siteCovs(x),
                         yearlySiteCovs=ysc, obsCovs=oc,
                         numPrimary=length(j))
})


# Select a subset of both sites i and observations j --------------------------
# i and j are numeric vectors
# many unmarkedFrames do not support this, I think


# i is as before and j is the obsNum to remove and corresponding y's
setMethod("[", c("unmarkedFrame","numeric", "numeric", "missing"),
		function(x, i, j)
{
    ## first remove sites
    umf <- x[i,]
    umf <- umf[,j]
    umf
})


# Select a subset of sites as a list-------------------------------------------

### list is a ragged array of indices (y's) to include for each site.
### Typically useful for multilevel boostrapping.
setMethod("[", c("unmarkedFrame","list", "missing", "missing"),
    function(x, i, j)
{
    m <- numSites(x)
    J <- R <- obsNum(x)
    o2y <- obsToY(x)
    if (!identical(o2y, diag(R)))
        stop("Ragged subsetting of unmarkedFrames is only valid for diagonal obsToY.")
    J <- ncol(o2y)
    if (m != length(i)) stop("list length must be same as number of sites.")
    siteCovs <- siteCovs(x)
    y <- cbind(.site=1:m, getY(x))
    obsCovs <- obsCovs(x)
    site_idx <- rep(1:m, each=R)
    stopifnot(length(site_idx) == nrow(obsCovs))

    oc <- lapply(1:m, function(ind){
      df <- obsCovs[site_idx==ind,,drop=FALSE]
      obs <- i[[ind]]
      if (length(obs) > R)
        stop("All elements of list must be less than or equal to R.")
      obs <- c(obs, rep(NA, R-length(obs)))
      df[obs,,drop=FALSE]
    })
    obsCovs <- do.call(rbind, oc)
    rownames(obsCovs) <- NULL

    y <- apply(y, 1, function(row) {
        site <- row[1]
        row <- row[-1]
        obs <- i[[site]]
        obs <- c(obs, rep(NA, R-length(obs)))
        row[obs]
        })

    if(!is.null(obsCovs(x))){
      obsCovs(x) <- obsCovs
    }
    x@y <- t(y)
    x
})


# Select subset of sites using logical vector----------------------------------

setMethod("[", c("unmarkedFrame", "logical", "missing", "missing"),
  function(x, i, j) {
  i <- which(i)
  x[i, ]
})

# Get first few rows of an unmarkedFrame---------------------------------------

setMethod("head", "unmarkedFrame", function(x, n) {
    if(missing(n)) n <- 10
    umf <- x[1:n,]
    umf
})

