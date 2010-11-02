
############# VALIDATION FUNCTIONS #############################################

validunmarkedFrame <- function(object) 
{
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
    if(J < 2) # matrices can have 0 columns
        errors <- c(errors, "y must have at least 2 columns")    
    if(!is.null(object@siteCovs))
        if(nrow(object@siteCovs) != M)
            errors <- c(errors, 
                "siteCovData does not have same size number of sites as y.")
    if(!is.null(obsCovs(object)) & !is.null(obsNum(object)))
        if(nrow(object@obsCovs) != M*obsNum(object))
            errors <- c(errors, "obsCovData does not have M*obsNum rows.")
    if(length(errors) == 0)
        TRUE
    else
        errors
}

############ DATA CLASSES ######################################################

# Class to hold data for analyses in unmarked.
setClass("unmarkedFrame",
    representation(y = "matrix",
        obsCovs = "optionalDataFrame",
        siteCovs = "optionalDataFrame",
        mapInfo = "optionalMapInfo",
        obsToY = "optionalMatrix"),
    validity = validunmarkedFrame)

## a class for multi-season data

setClass("unmarkedMultFrame",
    representation(numPrimary = "numeric",
        #data frame in site-major, year-minor order describing site-level covs
        yearlySiteCovs = "optionalDataFrame"),
    contains="unmarkedFrame")


## a class for distance sampling data
setClass("unmarkedFrameDS", 
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrame",
    validity = function(object) {
        errors <- character(0)
        J <- numY(object)
        db <- object@dist.breaks
        if(J != length(db) - 1)
            errors <- c(errors, "ncol(y) must equal length(dist.breaks)-1")
        if(db[1] != 0)
            errors <- c(errors, "dist.breaks[1] must equal 0")
        if(!is.null(obsCovs(object)))
            "obsCovs cannot be used with distsamp"
        if(length(errors) == 0) TRUE	
        else errors
        })


setClass("unmarkedFrameOccu",
		contains = "unmarkedFrame")


setClass("unmarkedFramePCount",
		contains = "unmarkedFrame")


setClass("unmarkedFrameMPois",
		representation(
			samplingMethod = "character",
			piFun = "character"),
		contains = "unmarkedFrame")


# for gmm aka gmultmix     
setClass("unmarkedFrameGMM", 
    representation(
        piFun = "character",
        samplingMethod = "character"),
    contains = "unmarkedMultFrame")    



################### CONSTRUCTORS ###############################################

# Constructor for unmarkedFrames.
unmarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo,
    obsToY) 
{
    if(class(obsCovs) == "list") {
        obsVars <- names(obsCovs)
        for(i in seq(length(obsVars))) {
            if(!(class(obsCovs[[i]]) %in% c("matrix", "data.frame")))
                stop("At least one element of obsCovs is not a matrix or data frame.")
            if(ncol(obsCovs[[i]]) != ncol(y) | nrow(obsCovs[[i]]) != nrow(y))
                stop("At least one matrix in obsCovs has incorrect number of dimensions.")
            }
        if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]]) #??
        obsCovs <- data.frame(lapply(obsCovs, function(x) as.vector(t(x))))
        }    
    if(("data.frame" %in% class(y)) | ("cast_matrix" %in% class(y))) 
        y <- as.matrix(y)    
    if(missing(obsToY)) obsToY <- NULL
    if(missing(mapInfo)) mapInfo <- NULL
	   
    umf <- new("unmarkedFrame", y = y, obsCovs = obsCovs, siteCovs = siteCovs, 
        mapInfo = mapInfo, obsToY = obsToY)
    return(umf)
}


unmarkedFrameDS <- function(y, siteCovs = NULL, dist.breaks, tlength, survey,
		unitsIn, mapInfo = NULL)
{
	if(missing(survey))
        stop("survey argument must be specified")
    if(missing(tlength) & survey == "point")
		tlength <- numeric(0)
	umfds <- new("unmarkedFrameDS", y = y, obsCovs = NULL,
			siteCovs = siteCovs, dist.breaks = dist.breaks, tlength = tlength,
			survey = survey, unitsIn = unitsIn,
			obsToY = matrix(1, 1, ncol(y)))
	return(umfds)
}



unmarkedFrameOccu <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo) 
{
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J), 
		mapInfo = mapInfo)
	umf <- as(umf, "unmarkedFrameOccu")
	umf
}

# This function constructs an unmarkedMultFrame object.
unmarkedMultFrame <- function(y, siteCovs = NULL, obsCovs = NULL, numPrimary,
	yearlySiteCovs = NULL) 
{
    J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J))
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
	
    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix", "data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != (ncol(y) / numPrimary) | 
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x) 
            as.vector(t(x))))
        }    

	umf@yearlySiteCovs <- yearlySiteCovs
	umf
}


# This function constructs an unmarkedMultFrame object.
unmarkedFrameGMM <- function(y, siteCovs = NULL, obsCovs = NULL, numPrimary,
	yearlySiteCovs = NULL, type, obsToY, piFun) 
{
    if(!missing(type)) {
      if(!type %in% c("removal", "double"))
        stop("if specifying type, it should either be 'removal' or 'double'")
      switch(type,
        removal = {
          obsToY <- matrix(1, ncol(y), ncol(y))
				  obsToY[col(obsToY) < row(obsToY)] <- 0
				  piFun <- "removalPiFun"
          },
        double = {
				  obsToY <- matrix(c(1, 0, 0, 1, 1, 1), 2, 3)
				  piFun <- "doublePiFun"
          })
      }
    else {
        type <- "userDefined" 
        if(missing(obsToY)) 
            stop("obsToY is required for gmultmix data with no specified type.")
		}
       
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    J <- ncol(y) / numPrimary
    if(class(yearlySiteCovs) == "list") {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!(class(yearlySiteCovs[[i]]) %in% c("matrix", "data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != numPrimary | 
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        if(is.null(obsNum)) obsNum <- ncol(obsCovs[[1]])
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x) 
            as.vector(t(x))))
        }    
    umf@yearlySiteCovs <- yearlySiteCovs
    umf <- as(umf, "unmarkedFrameGMM")
    umf@piFun <- piFun
    umf@samplingMethod <- type
    umf
}



unmarkedFramePCount <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo) 
{
	J <- ncol(y)
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J), 
		mapInfo = mapInfo)
	umf <- as(umf, "unmarkedFramePCount")
	umf
}



unmarkedFrameMPois <- function(y, siteCovs = NULL, obsCovs = NULL, type, obsToY, 
	mapInfo, piFun) 
{
	if(!missing(type)) {
		switch(type,
			removal = {
				obsToY <- matrix(1, ncol(y), ncol(y))
				obsToY[col(obsToY) < row(obsToY)] <- 0
				piFun <- "removalPiFun"
			},
			double = {
				obsToY <- matrix(c(1, 0, 0, 1, 1, 1), 2, 3)
				piFun <- "doublePiFun"
			})
		}
	else {
		if(missing(obsToY)) 
			stop("obsToY is required for multinomial-Poisson data with no specified type.")
		type <- "userDefined"
		}
	umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY, 
		mapInfo = mapInfo)
	umf <- as(umf, "unmarkedFrameMPois")
	umf@piFun <- piFun
	umf@samplingMethod <- type
	umf
}

################ SHOW METHODS ##################################################


setMethod("show", "unmarkedFrame",
        function(object) {
            df <- as(object, "data.frame")
            cat("Data frame representation of unmarkedFrame object.\n")
            print(df)
        })


setMethod("show", "unmarkedMultFrame",
        function(object) {
            df <- as(object, "data.frame")
            ysc <- yearlySiteCovs(object)
            if(is.null(ysc)) { 
                cat("Data frame representation of unmarkedFrame object.\n")
                print(df)
                }
            else {
                T <- object@numPrimary
                yscwide <- data.frame(lapply(ysc, matrix, ncol=T))
                df <- data.frame(df, yscwide)
                cat("Data frame representation of unmarkedFrame object.\n")
                print(df)
                }
        })


############################ EXTRACTORS ########################################

# Extractor for site level covariates
setGeneric("siteCovs", function(object,...) standardGeneric("siteCovs"))
setMethod("siteCovs", "unmarkedFrame",
		function(object) {
			return(object@siteCovs)
		})

setGeneric("yearlySiteCovs", 
	function(object,...) standardGeneric("yearlySiteCovs"))
setMethod("yearlySiteCovs", "unmarkedMultFrame",
		function(object) {
			return(object@yearlySiteCovs)
		})

setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))
setMethod("obsCovs", "unmarkedFrame", 
		function(object, matrices = FALSE) {
			M <- numSites(object)
			R <- obsNum(object)
			if(matrices) {
				value <- list()
				for(i in seq(length=length(object@obsCovs))){
					value[[i]] <- matrix(object@obsCovs[,i], M, R, byrow = TRUE)
				}
				names(value) <- names(object@obsCovs)
			} else {
				value <- object@obsCovs
			}
			return(value)
		})


setGeneric("obsNum", function(object) standardGeneric("obsNum"))
setMethod("obsNum", "unmarkedFrame", function(object) nrow(object@obsToY))


setGeneric("numSites", function(object) standardGeneric("numSites"))
setMethod("numSites", "unmarkedFrame", function(object) nrow(object@y))


setGeneric("numY", function(object) standardGeneric("numY"))
setMethod("numY", "unmarkedFrame", function(object) ncol(object@y))


setGeneric("obsToY", function(object) standardGeneric("obsToY"))
setMethod("obsToY", "unmarkedFrame", function(object) object@obsToY)


setGeneric("obsCovs<-", function(object, value) standardGeneric("obsCovs<-"))
setReplaceMethod("obsCovs", "unmarkedFrame", function(object, value) {
			if(identical(class(object)[1], "unmarkedFrameDS"))
			     stop("unmarkedFrameDS objects cannot have obsCovs")     
            object@obsCovs <- as.data.frame(value)
			object
		})


setGeneric("siteCovs<-", function(object, value) standardGeneric("siteCovs<-"))
setReplaceMethod("siteCovs", "unmarkedFrame", function(object, value) {
			object@siteCovs <- as.data.frame(value)
			object
		})


setGeneric("yearlySiteCovs<-", 
	function(object, value) standardGeneric("yearlySiteCovs<-"))
setReplaceMethod("yearlySiteCovs", "unmarkedMultFrame", function(object, value) {
			object@yearlySiteCovs <- as.data.frame(value)
			object
		})


setGeneric("obsToY<-", function(object, value) standardGeneric("obsToY<-"))
setReplaceMethod("obsToY", "unmarkedFrame", function(object, value) {
			object@obsToY <- value
			object
		})



setGeneric("getY", function(object) standardGeneric("getY"))
setMethod("getY", "unmarkedFrame", function(object) object@y)


setGeneric("coordinates", function(object) standardGeneric("coordinates"))
setMethod("coordinates", "unmarkedFrame",
		function(object) {
			object@mapInfo@coordinates
		})


setGeneric("projection", function(object) standardGeneric("projection"))
setMethod("projection", "unmarkedFrame",
		function(object) {
			object@mapInfo@projection
		})

################################### SUMMARY METHODS ############################


setMethod("summary", "unmarkedFrame",
	function(object,...) {
		cat("unmarkedFrame Object\n\n")
		cat(nrow(object@y), "sites\n")
		cat("Maximum number of observations per site:",obsNum(object),"\n")
			mean.obs <- mean(rowSums(!is.na(getY(object))))
		cat("Mean number of observations per site:",round(mean.obs,2),"\n")
		cat("Sites with at least one detection:", 
			sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))), 
				"\n\n")
		cat("Tabulation of y observations:")
		print(table(object@y, exclude=NULL))
		if(!is.null(object@siteCovs)) {
			cat("\nSite-level covariates:\n")
			print(summary(object@siteCovs))
			}
		if(!is.null(object@obsCovs)) {
			cat("\nObservation-level covariates:\n")
			print(summary(object@obsCovs))
			}
	})

	

setMethod("summary", "unmarkedFrameDS", 
	function(object, ...) 
{
	cat("unmarkedFrameDS Object\n\n")
	cat(object@survey, "-transect survey design", "\n", sep="")
	cat(paste("Distance class cutpoints (", object@unitsIn, "): ", sep=""), 
		object@dist.breaks, "\n\n")
	cat(nrow(object@y), "sites\n")
	cat("Maximum number of distance classes per site:", ncol(getY(object)), "\n")
		mean.dc <- mean(rowSums(!is.na(getY(object))))
	cat("Mean number of distance classes per site:", round(mean.dc, 2), "\n")
	cat("Sites with at least one detection:", 
		sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))), "\n\n")
	cat("Tabulation of y observations:")
	print(table(object@y, exclude=NULL))
	if(!is.null(object@siteCovs)) {
		cat("\nSite-level covariates:\n")
		print(summary(object@siteCovs))
		}
	if(!is.null(object@obsCovs)) {
		warning("Observation-level covariates cannot be used by distsamp()")
		}
})




setMethod("summary", "unmarkedMultFrame",
	function(object,...) {
		cat("unmarkedFrame Object\n\n")
		cat(nrow(object@y), "sites\n")
		cat("Maximum number of observations per site:",obsNum(object),"\n")
			mean.obs <- mean(rowSums(!is.na(getY(object))))
		cat("Mean number of observations per site:",round(mean.obs,2),"\n")
		cat("Number of primary survey periods:", object@numPrimary, "\n")
		cat("Number of secondary survey periods:", 
            obsNum(object) / object@numPrimary, "\n")
        cat("Sites with at least one detection:", 
			sum(apply(getY(object), 1, function(x) any(x > 0, na.rm=TRUE))), 
				"\n\n")
		cat("Tabulation of y observations:")
		print(table(object@y, exclude=NULL))
		if(!is.null(object@siteCovs)) {
			cat("\nSite-level covariates:\n")
			print(summary(object@siteCovs))
			}
		if(!is.null(object@obsCovs)) {
			cat("\nObservation-level covariates:\n")
			print(summary(object@obsCovs))
			}
		if(!is.null(object@yearlySiteCovs)) {
			cat("\nYearly-site-level covariates:\n")
			print(summary(object@yearlySiteCovs))
			}
	})




################################# PLOT METHODS #################################
# TODO:  come up with nice show/summary/plot methods for each of these data types.

setMethod("plot", c(x="unmarkedFrame", y="missing"),
	function (x, y, panels = 1, colorkey, strip=FALSE, 
    ylab="Site", xlab="Observation", ...)
{
    y <- getY(x)
    ym <- max(y, na.rm=TRUE)
    M <- nrow(y)
    J <- ncol(y)
    y <- as.data.frame(y)
    colnames(y) <- paste("obs",1:J)
    y$site <- 1:M
    sites.per.panel <- M/panels
    y$group <- as.factor(round(seq(1,panels,length=M)))
    y2 <- melt(y, #measure.var = c("V1", "V2", "V3"),
        id.var=c("site","group"))
    if(missing(colorkey))
        colorkey <- list(at=0:(ym+1), labels=list(labels=as.character(0:ym), 
            at=(0:ym)+0.5))    
    levelplot(value ~ variable*site | group, y2,
        scales=list(relation="free", x=list(labels=1:J)), 
        colorkey=colorkey, strip=strip, xlab=xlab, ylab=ylab, ...)
})


setMethod("hist", "unmarkedFrameDS", function(x, ...)
{
    y <- getY(x)
    dbreaks <- x@dist.breaks
    nb <- length(dbreaks)
    mids <- (dbreaks[-1] - dbreaks[-nb]) / 2 + dbreaks[-nb]
        distances <- rep(mids, times=colSums(y))
    hist(distances, breaks=dbreaks, ...)
})



#setMethod("plot", c(x="unmarkedFrameOccu", y="missing"),
#		function(x) {
#			if(is.null(x@mapInfo)) stop("mapInfo is required to plot an unmarkedFrameOccu object.")
#			y <- getY(x)
#			## get sites w/ at least one pos
#			y <- as.factor(rowSums(y, na.rm = TRUE) > 0)
#			levels(y) <- c("non-detection", "detection")
#			siteCovs <- siteCovs(x)
#			coords <- coordinates(x)
#			if(is.null(x@mapInfo@projection)) {
#				proj <- list(x = coords[,1], y = coords[,2])	
#			} else {
#				proj <- mapproject(x = coords[,1], y = coords[,2], projection = x@mapInfo@projection,
#						parameters = x@mapInfo@parameters, orientation = x@mapInfo@orientation)
#			}
#			p <- qplot(x = proj$x, y = proj$y, colour = y, xlab = "longitude", ylab = "latitude")
#			if(!is.null(x@mapInfo@projection)) {
#				p + coord_map(project = x@mapInfo@projection)
#			} else {
#				p
#			}
#		})
#

#setMethod("plot", c(x="unmarkedFramePCount", y="missing"),
#		function(x) {
#			if(is.null(x@mapInfo)) stop("mapInfo is required to plot an unmarkedFramePCount object.")
#			y <- getY(x)
#			## plot maximums
#			y <- apply(y, 1, max, na.rm = TRUE)
#			siteCovs <- siteCovs(x)
#			coords <- coordinates(x)
#			if(is.null(x@mapInfo@projection)) {
#				proj <- list(x = coords[,1], y = coords[,2])	
#			} else {
#				proj <- mapproject(x = coords[,1], y = coords[,2], projection = x@mapInfo@projection,
#					parameters = x@mapInfo@parameters, orientation = x@mapInfo@orientation)
#			}
#			p <- qplot(x = proj$x, y = proj$y, colour = y, xlab = "longitude", ylab = "latitude")
#			if(!is.null(x@mapInfo@projection)) {
#				p + coord_map(project = x@mapInfo@projection)
#			} else {
#				p
#			}
#		})
################################# SELECTORS ####################################

# i is the vector of sites to extract

setMethod("[", c("unmarkedFrame", "numeric", "missing", "missing"),
	function(x, i) {  
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
            obsCovs <- cbind(.site=rep(1:M, each = R), obsCovs(x))
            obsCovs <- ldply(i, function(site) {
              subset(obsCovs, .site == site)
            })
            obsCovs$.site <- NULL
          }
          umf <- x
          umf@y <- y
          umf@siteCovs <- siteCovs
          umf@obsCovs <- obsCovs
          umf
	})

## remove obs only
setMethod("[", c("unmarkedFrame", "missing", "numeric", "missing"),
		function(x, i, j) {  
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

# i is as before and j is the obsNum to remove and corresponding y's
setMethod("[", c("unmarkedFrame","numeric", "numeric", "missing"),
		function(x, i, j) {  
			## first remove sites
			umf <- x[i,]
			umf <- umf[,j]
			umf
		})


### list is a ragged array of indices (y's) to include for each site.
### Typically useful for multilevel boostrapping.
setMethod("[", c("unmarkedFrame","list", "missing", "missing"),
function(x, i) {  
  m <- numSites(x)
  J <- R <- obsNum(x)
  o2y <- obsToY(x)
  if (!identical(o2y, diag(R))) stop("Ragged subsetting of unmarkedFrames is only valid for diagonal obsToY.")
  J <- ncol(o2y)
  if (m != length(i)) stop("list length must be same as number of sites.")
  siteCovs <- siteCovs(x)
  y <- cbind(.site=1:m, getY(x))
  obsCovs <- cbind(.site=rep(1:m, each=R), obsCovs(x))

  obsCovs <- ddply(obsCovs, ~.site, function(df) {
    site <- df$.site[1]
    obs <- i[[site]]
    if (length(obs) > R) stop("All elements of list must be less than or equal to R.")
    obs <- c(obs, rep(NA, R-length(obs)))
    df[obs,]
  })
  obsCovs$.site <- NULL

  y <- apply(y, 1, function(row) {
    site <- row[1]
    row <- row[-1]
    obs <- i[[site]]
    obs <- c(obs, rep(NA, R-length(obs)))
    row[obs]
  })

  obsCovs(x) <- obsCovs
  x@y <- t(y)
  x
})




## for multframes, must remove years at a time
setMethod("[", c("unmarkedMultFrame", "missing", "numeric", "missing"),
		function(x, i, j) {  
			J <- obsNum(x)/x@numPrimary
			obs <- rep(1:x@numPrimary, each = J)
			numPrimary <- length(j)
			j <- which(!is.na(match(obs,j)))
			u <- callNextMethod(x, i, j)
			u@numPrimary <- numPrimary
			u
		})


setMethod("head", "unmarkedFrame",
		function(x, n) {
			if(missing(n)) n <- 10
			umf <- x[1:n,]
			umf
		})
		
############################### COERCION #######################################

setAs("data.frame", "unmarkedFrame", function(from) {
			umf <- formatWide(from)
			umf
		})

setAs("unmarkedFrame", "data.frame", function(from) {
			obsCovs <- obsCovs(from)
			siteCovs <- siteCovs(from)
			y <- getY(from)
			colnames(y) <- paste("y",1:ncol(y),sep=".")
			if(is.null(obsToY(from))) {
				obsNum <- ncol(y)
			} else {
				obsNum <- obsNum(from)
			}
			if(is.null(siteCovs)) siteCovs <- matrix(0,nrow(y),0)
			if(is.null(obsCovs)) {
				obsCovs <- matrix(0,nrow(y),0)
			} else {
				obsCovs <- data.frame(lapply(obsCovs, 
					function(x) matrix(x, nrow(y), obsNum,byrow=T)))
			}
			df <- data.frame(y, siteCovs, obsCovs)
			df
		})
		
		
		

			
		
