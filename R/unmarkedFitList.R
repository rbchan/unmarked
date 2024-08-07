#Get y (corrected for missing values) based on data and formulas
#Approach for some fit functions is different, so this is now a method
setGeneric("fl_getY", function(fit, ...) standardGeneric("fl_getY"))

setMethod("fl_getY", "unmarkedFit", function(fit, ...){
  getDesign(getData(fit), fit@formula)$y
})

setMethod("fl_getY", "unmarkedFitOccuMulti", function(fit, ...){
  maxOrder <- fit@call$maxOrder
  if(is.null(maxOrder)) maxOrder <- length(fit@data@ylist)
  getDesign(getData(fit), fit@detformulas, fit@stateformulas, maxOrder)$y
})

setMethod("fl_getY", "unmarkedFitOccuMS", function(fit, ...){
  getDesign(getData(fit), fit@psiformulas, fit@phiformulas,
            fit@detformulas, fit@parameterization)$y
})

setMethod("fl_getY", "unmarkedFitOccuFP", function(fit, ...){
  getDesign(getData(fit), fit@detformula, fit@FPformula,
            fit@Bformula, fit@stateformula)$y
})

setClass("unmarkedFitList",
    representation(fits = "list"),
    validity = function(object) {
        fl <- object@fits
        umf1 <- getData(fl[[1]])
        y1 <- fl_getY(fl[[1]])
        dataTest <- sapply(fl, function(x) isTRUE(all.equal(umf1, getData(x))))
        yTest <- sapply(fl, function(x) isTRUE(all.equal(y1, fl_getY(x))))
        if(!all(dataTest)) {
            stop("Data are not the same among models. Make sure you use the same unmarkedFrame object for all models.")
            }
        else if(!all(yTest)) {
            stop("Data are not the same among models due to missing covariate values. Consider removing NAs before analysis.")
            }
        TRUE
        }
    )


# constructor of unmarkedFitList objects
fitList <- function(..., fits, autoNames=c("object","formula")) {
  autoNames <- match.arg(autoNames)
  if(length(list(...)) > 0 & !missing(fits))
    stop("Do not use both the '...' and 'fits' arguments")
  if(missing(fits)) {
    fits <- list(...)
    isList <- sapply(fits, function(x) is.list(x))
    if(sum(isList) > 1)
      stop("Specify models as common-seperated objects, or use fits = 'mylist'")
    if(isList[1L]) {
      warning("If supplying a list of fits, use fits = 'mylist'")
      fits <- fits[[1L]] 	# This is allowed for back-compatability.
    }
    if(is.null(names(fits))) {
      message("Your list was unnamed, so model names were added automatically")
      if(autoNames=="formula"){
        if(inherits(fits[[1]], c("unmarkedFitOccuMulti", "unmarkedFitOccuMS"))){
          warning("simple formula not available, naming based on object instead", call.=FALSE)
          autoNames <- "object"
        } else {
          names(fits) <- lapply(fits, function(x) gsub(" ", "", as.character(deparse(x@formula))))
        }
      }
      if(autoNames=="object"){
        c <- match.call(expand.dots = FALSE)
        names(fits) <- as.character(c[[2]])
      }
    }
  }
  if(is.null(names(fits))) {
    message("Your list was unnamed, so model names were added automatically")
    if(autoNames=="formula"){
        if(inherits(fits[[1]], c("unmarkedFitOccuMulti", "unmarkedFitOccuMS"))){
        warning("simple formula not available, naming as numbers instead", call.=FALSE)
        autoNames <- "object"
      } else {
        names(fits) <- lapply(fits, function(x) gsub(" ", "", as.character(deparse(x@formula))))
      }
    }
    if(autoNames=="object"){
      names(fits) <- as.character(1:(length(fits)))
    }
  }
  umfl <- new("unmarkedFitList", fits=fits)
  return(umfl)
}

setMethod("summary", "unmarkedFitList", function(object) {
    fits <- object@fits
    for(i in 1:length(fits))
        summary(fits[[i]])
    })


setMethod("coef", "unmarkedFitList", function(object)
{
    fits <- object@fits
    coef.list <- lapply(fits, coef)
    coef.names <- unique(unlist(lapply(coef.list, names)))
    coef.out <- matrix(NA, length(fits), length(coef.names))
    colnames(coef.out) <- coef.names
    if(!is.null(names(fits)))
        rownames(coef.out) <- names(fits)
    for(i in 1:length(coef.list))
        coef.out[i, names(coef.list[[i]])] <- coef.list[[i]]
    return(coef.out)
})


setMethod("SE", "unmarkedFitList", function(obj)
{
    fits <- obj@fits
    se.list <- lapply(fits, function(x) {
        tr <- try(SE(x))
        if(class(tr)[1] == "try-error")
            return(rep(NA, length(coef(x))))
        else
            return(tr)
        })
    se.names <- unique(unlist(lapply(se.list, names)))
    se.out <- matrix(NA, length(fits), length(se.names))
    colnames(se.out) <- se.names
    if(!is.null(names(fits)))
        rownames(se.out) <- names(fits)
    for(i in 1:length(se.list))
        se.out[i, names(se.list[[i]])] <- se.list[[i]]
    return(se.out)
})




setMethod("predict", "unmarkedFitList", function(object, type, newdata=NULL,
    backTransform = TRUE, appendData = FALSE, level=0.95) {
        fitList <- object@fits
        ese <- lapply(fitList, predict, type = type, newdata = newdata,
            backTransform = backTransform, level=level)

        if(inherits(newdata, "Raster")){
          if(!requireNamespace("raster")) stop("raster package is required")
          ese <- lapply(ese, raster::as.matrix)
        }

        E <- sapply(ese, function(x) x[,"Predicted"])
        SE <- sapply(ese, function(x) x[,"SE"])
        lower <- sapply(ese, function(x) x[,"lower"])
        upper <- sapply(ese, function(x) x[,"upper"])
        ic <- sapply(fitList, slot, "AIC")
        deltaic <- ic - min(ic)
        wts <- exp(-deltaic / 2)
        wts <- wts / sum(wts)
        parav <- as.numeric(E %*% wts)
        seav <- as.numeric(sqrt(SE^2 + (E - parav)^2) %*% wts)
        out <- data.frame(Predicted = parav, SE = seav)
        out$lower <- as.numeric(lower %*% wts)
        out$upper <- as.numeric(upper %*% wts)

        if(inherits(newdata, "Raster")){
          E.mat <- matrix(out[,1], dim(newdata)[1], dim(newdata)[2], byrow=TRUE)
          E.raster <- raster::raster(E.mat)
          raster::extent(E.raster) <- raster::extent(newdata)
          out.rasters <- list(E.raster)
          for(i in 2:ncol(out)) {
            i.mat <- matrix(out[,i], dim(newdata)[1], dim(newdata)[2], byrow=TRUE)
            i.raster <- raster::raster(i.mat)
            raster::extent(i.raster) <- raster::extent(newdata)
            out.rasters[[i]] <- i.raster
          }
          out.stack <- raster::stack(out.rasters)
          names(out.stack) <- colnames(out)
          raster::crs(out.stack) <- raster::crs(newdata)
          return(out.stack)
        }

        if(appendData) {
            if(missing(newdata))
                newdata <- getData(object@fits[[1]])
            out <- data.frame(out, newdata)
            }
        return(out)
    })






# Condition number
cn <- function(object) {
    h <- hessian(object)
    if(is.null(h)) return(NA)
    if(any(is.na(h))) return(NA)
        else {
   	        ev <- eigen(h)$value
   	        return(max(ev) / min(ev))
   	        }
   	}



# R-squared index from Nagelkerke (1991)
nagR2 <- function(fit, nullfit)
{
    n <- sampleSize(fit)
    devI <- 2 * fit@negLogLike
    devN <- 2 * nullfit@negLogLike
    r2 <- 1 - exp((devI - devN) / n)
    r2max <- 1 - exp(-1 * devN / n)
    return(r2 / r2max)
}



setGeneric("modSel",
        def = function(object, ...) {
            standardGeneric("modSel")
            }
        )

setClass("unmarkedModSel",
    representation(
        Full = "data.frame",
        Names = "matrix"
        )
    )



# Model selection results from an unmarkedFitList
setMethod("modSel", "unmarkedFitList",
	function(object, nullmod=NULL)
{
    if (!is.character(nullmod) && !is.null(nullmod)) {
        stop("nullmod must be character name of null model fit in the fitlist.")
        }
    fits <- object@fits
    estList <- lapply(fits, coef, altNames=TRUE)
    seList <- lapply(fits, function(x) {
		se <- tryCatch(sqrt(diag(vcov(x, altNames=TRUE))),
			error=function(e) simpleError("Hessian is singular."))
        if(identical(class(se)[1], "simpleError")) {
            cat(se$message, fill=TRUE)
            se <- rep(NA, length(coef(x)))
            }
        return(se)
        })
    eNames <- sort(unique(unlist(sapply(estList, names))))
    seNames <- paste("SE", eNames, sep="")
    eseNames <- character(l <- length(c(eNames, seNames)))
    eseNames[seq(1, l, by=2)] <- eNames
    eseNames[seq(2, l, by=2)] <- seNames
    cNames <- c("model", "formula", eseNames)
    out <- data.frame(matrix(NA, ncol=length(cNames), nrow=length(fits)))
    colnames(out) <- cNames
    out$model <- names(fits)
    out$formula <- sapply(fits, function(x) {
          f <- as.character(x@formula)
          f <- paste(f[2], "~", f[3])
          f
        })
    for(i in 1:length(eNames)) {
        out[,eNames[i]] <- sapply(estList, function(x) x[eNames[i]])
        out[,seNames[i]] <- sapply(seList, function(x) x[eNames[i]])
        }
    out$Converge <- sapply(fits, function(x) x@opt$convergence)
    out$CondNum <- sapply(fits, function(x) cn(x))
    out$negLogLike <- sapply(fits, function(x) x@negLogLike)
    out$nPars <- sapply(fits, function(x) length(coef(x)))
    out$n <- sapply(fits, function(x) sampleSize(x))
    out$AIC <- sapply(fits, function(x) x@AIC)
    out$delta <- out$AIC - min(out$AIC)
    out$AICwt <- exp(-out$delta / 2)
    out$AICwt <- out$AICwt / sum(out$AICwt)
    out$Rsq <- NA
    if(!is.null(nullmod)) {
          if (is.na(match(nullmod, names(fits)))) {
            stop(paste("No fit named", nullmod, "was found in fits."))
          }
          nullmod <- fits[[nullmod]]
          out$Rsq <- sapply(fits, nagR2, nullmod)
        }
    out <- out[order(out$AIC),]
    out$cumltvWt <- cumsum(out$AICwt)
    msout <- new("unmarkedModSel", Full = out,
        Names = rbind(Coefs = eNames, SEs = seNames))
    return(msout)
})




setAs("unmarkedModSel", "data.frame", function(from) {
    out <- from@Full
    out
})



setMethod("show", "unmarkedModSel", function(object)
{
    out <- as(object, "data.frame")
    rownames(out) <- out$model
    out <- out[,c('nPars', 'AIC', 'delta', 'AICwt', 'cumltvWt', 'Rsq')]
    if (all(is.na(out$Rsq))) out$Rsq <- NULL
    print(format(out, digits=2, nsmall=2))
})



setMethod("coef", "unmarkedModSel", function(object)
{
    coefNames <- object@Names["Coefs",]
    msdf <- as(object, "data.frame")
    rownames(msdf) <- msdf$model
    out <- msdf[,coefNames]
    out
})


setMethod("SE", "unmarkedModSel", function(obj)
{
    seNames <- obj@Names["SEs",]
    msdf <- as(obj, "data.frame")
    rownames(msdf) <- msdf$model
    out <- msdf[,seNames]
    out
})


