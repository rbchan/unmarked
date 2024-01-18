
# ------------------------ VALIDATION FUNCTIONS --------------------------

validunmarkedFrame <- function(object) {
    errors <- character(0)
    M <- nrow(object@y)
    J <- ncol(object@y)
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

# --------------------------- DATA CLASSES -------------------------------

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
        #data frame in site-major, year-minor order describing siteCovs
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

setClass("unmarkedFrameOccuFP",
         representation(
           type = "numeric"),
         contains = "unmarkedFrame")

#Multi-state occupancy
setClass('unmarkedFrameOccuMS',
         representation(
            numStates = "numeric",
            phiOrder = "list"),
         contains = "unmarkedMultFrame")

#Time-to-detection occupancy
setClass("unmarkedFrameOccuTTD",
         representation(
            surveyLength = "matrix"),
          contains = "unmarkedMultFrame")

setClass("unmarkedFrameOccuMulti",
    representation(ylist = "list", fDesign = "matrix"),
    contains = "unmarkedFrame",
    validity = function(object) {
      errors <- character(0)
      M <- nrow(object@y)
      J <- ncol(object@y)
      Ms <- sapply(object@ylist,nrow)
      Js <- sapply(object@ylist,ncol)
      if(length(unique(Ms)) != 1)
        errors <- c(errors, "All species must have same number of sites")
      if(length(unique(Js)) != 1)
        errors <- c(errors, "All species must have same number of observations")
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
      })

setClass("unmarkedFramePCount",
		contains = "unmarkedFrame")


setClass("unmarkedFrameMPois",
		representation(
			samplingMethod = "character",
			piFun = "character"),
		contains = "unmarkedFrame")


setClass("unmarkedFrameG3",
         contains = "unmarkedMultFrame")


setClass("unmarkedFramePCO",
         representation(primaryPeriod = "matrix"),
         contains = "unmarkedMultFrame")

setClass("unmarkedFrameGMM",
    representation(
        piFun = "character",
        samplingMethod = "character"),
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGDS",
    representation(
        dist.breaks = "numeric",
        tlength = "numeric",
        survey = "character",
        unitsIn = "character"),
    contains = "unmarkedFrameG3")

setClass("unmarkedFrameGPC",
    contains = "unmarkedFrameG3")


setClass("unmarkedFrameMMO",
         representation(primaryPeriod = "matrix"),
         contains = "unmarkedFrameGMM")

## Andy 12/27/2015
setClass("unmarkedFrameDSO",
         representation(primaryPeriod = "matrix"),
         contains = "unmarkedFrameGDS")


# ------------------------------- CONSTRUCTORS ---------------------------

#Convert covs provided as list of matrices/dfs to data frame
covsToDF <- function(covs, name, obsNum, numSites){
  if(is.null(covs)) return(covs)
  if(inherits(covs, "data.frame")){
    if(nrow(covs) != (obsNum * numSites)){
      stop("Incorrect number of rows in ", name, " data frame", call.=FALSE)
    }
    return(covs)
  }
  
  if(is.null(names(covs)) | any(is.na(names(covs))) | any(names(covs)=="")){
    stop("All elements of list provided to ", name, " argument must be named", call.=FALSE)
  }
  lapply(covs, function(x){
    if(!inherits(x, c("matrix", "data.frame")))
      stop(paste("At least one element of", name, "is not a matrix or data frame."))
    if(ncol(x) != obsNum | nrow(x) != numSites)
      stop(paste("At least one element of", name, "has incorrect number of dimensions."))
  })

  data.frame(lapply(covs, function(x) as.vector(t(x))))
}

# Constructor for unmarkedFrames.
unmarkedFrame <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo,
                          obsToY) {
    if(!missing(obsToY)){
        obsNum <- nrow(obsToY)
    } else {
        obsNum <- ncol(y)
    }

    if(is.null(obsNum) & inherits(obsCovs, "list"))
      obsNum <- ncol(obsCovs[[1]]) #??

    obsCovs <- covsToDF(obsCovs, "obsCovs", obsNum, nrow(y))

    if(inherits(y, c("data.frame", "cast_matrix")))
        y <- as.matrix(y)
    if(missing(obsToY)) obsToY <- NULL
    if(missing(mapInfo)) mapInfo <- NULL

    umf <- new("unmarkedFrame", y = y, obsCovs = obsCovs, siteCovs = siteCovs,
        mapInfo = mapInfo, obsToY = obsToY)
    umf <- umf_to_factor(umf)
    return(umf)
}


unmarkedFrameDS <- function(y, siteCovs = NULL, dist.breaks, tlength,
                            survey, unitsIn, mapInfo = NULL)
{
    if(missing(survey))
        stop("survey argument must be specified")
    if(missing(tlength) & survey == "point")
        tlength <- rep(NA_real_, nrow(y))
    if((survey=="line") & (length(tlength) != nrow(y)))
        stop("tlength should be a vector with length(tlength)==nrow(y)")
    if(length(dist.breaks) != ncol(y)+1){
      stop(paste("dist.breaks should have",ncol(y)+1,"values"))
    }
    umfds <- new("unmarkedFrameDS", y = y, obsCovs = NULL,
                 siteCovs = siteCovs, dist.breaks = dist.breaks,
                 tlength = tlength, survey = survey, unitsIn = unitsIn,
                 obsToY = matrix(1, 1, ncol(y)))
    umfds <- umf_to_factor(umfds)
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

unmarkedFrameOccuFP <- function(y, siteCovs = NULL, obsCovs = NULL, type, mapInfo)
{
  J <- ncol(y)
  umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J),
                       mapInfo = mapInfo)
  umf <- as(umf, "unmarkedFrameOccuFP")
  umf@type <- type
  umf
}

unmarkedFrameOccuMulti <- function(y, siteCovs = NULL, obsCovs = NULL,
                                   maxOrder, mapInfo = NULL)
{
  ylist <- y

  #Fix non-matching missing observations in ylist
  na_list <- lapply(ylist, function(x) is.na(x))
  sum_na <- sapply(na_list, sum)

  if(any(sum_na[1] != sum_na)){
    warning("Missing obs do not match among species. Converting some obs to NAs so they match.")
  }

  any_na <- Reduce(function(x,y) x|y, na_list)
  ylist <- lapply(ylist, function(x){
                      x[any_na] <- NA
                      x})

  y <- ylist[[1]]
  J <- ncol(y)
  if(is.null(names(ylist)))
    names(ylist) <- paste('sp',1:length(ylist),sep='')

  obsCovs <- covsToDF(obsCovs, "obsCovs", J, nrow(y))

  #f design matrix guide
  S <- length(ylist)
  z <- expand.grid(rep(list(1:0),S))[,S:1]
  colnames(z) <- names(ylist)

  if(missing(maxOrder)) maxOrder <- S
  if(maxOrder == 1){
    fDesign <- as.matrix(z)
  } else {
    fDesign <- model.matrix(as.formula(paste0("~.^",maxOrder,"-1")),z)
  }

  attr(fDesign,'assign') <- NULL
  zinds <- apply(z,1,function(x) paste(x,collapse=''))
  rownames(fDesign) <- paste('psi[',zinds,']',sep='')
  colnames(fDesign) <- paste('f',1:ncol(fDesign),'[',
                             colnames(fDesign),']',sep='')

  umfmo <- new("unmarkedFrameOccuMulti", y=y, ylist = ylist, fDesign=fDesign,
               obsCovs = obsCovs, siteCovs = siteCovs, obsToY = diag(J))
  umfmo <- umf_to_factor(umfmo)
  return(umfmo)
}


unmarkedFramePCount <- function(y, siteCovs = NULL, obsCovs = NULL, mapInfo)
{
    J <- ncol(y)
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J),
        mapInfo = mapInfo)
    umf <- as(umf, "unmarkedFramePCount")
    umf
}



unmarkedFrameMPois <- function(y, siteCovs = NULL, obsCovs = NULL, type,
    obsToY, mapInfo, piFun)
{
    if(!missing(type)) {
        switch(type,
            removal = {
                obsToY <- matrix(1, ncol(y), ncol(y))
                obsToY[col(obsToY) < row(obsToY)] <- 0
                #obsToY <- diag(ncol(y))
                #obsToY[upper.tri(obsToY)] <- 1
                piFun <- "removalPiFun"
                },
            double = {
                #obsToY <- matrix(c(1, 0, 0, 1, 1, 1), 2, 3)
                obsToY <- matrix(1, 2, 3)
                piFun <- "doublePiFun"
                },
            depDouble = {
              obsToY <- matrix(1, 2, 2)
              piFun <- "depDoublePiFun"
              if(ncol(y) != 2){
                stop("y must have exactly 2 columns when type = 'depDouble'",
                     call.=FALSE)
              }
              })
    } else {
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



# This function constructs an unmarkedMultFrame object.
unmarkedMultFrame <- function(y, siteCovs = NULL, obsCovs = NULL,
                              numPrimary, yearlySiteCovs = NULL)
{
    J <- ncol(y)
	  umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J))
    umf <- as(umf, "unmarkedMultFrame")
    if(ncol(umf@y) %% numPrimary != 0){
      stop("Unequal number of secondary periods in each primary period", call.=FALSE)
    }
    umf@numPrimary <- numPrimary
    umf@yearlySiteCovs <- covsToDF(yearlySiteCovs, "yearlySiteCovs",
                                   numPrimary, nrow(y))
    umf <- umf_to_factor(umf)
    umf
}


# Construct a unmarkedFrameOccuMS object (multistate occupancy)
unmarkedFrameOccuMS <- function(y, siteCovs = NULL, obsCovs = NULL,
                                numPrimary = 1, yearlySiteCovs = NULL)
{

  umf <- unmarkedMultFrame(y, siteCovs, obsCovs, numPrimary, yearlySiteCovs)
  umf <- as(umf, "unmarkedFrameOccuMS")
  umf@numStates <- max(y,na.rm=T) + 1
  if(umf@numStates<3){
    stop("<3 occupancy states detected. Use occu() or colext() instead.")
  }

  #Create guide for phi formula order
  S <- umf@numStates

  #Multinomial
  vals <- paste0('phi[',rep(0:(S-1),each=S),'->',rep(0:(S-1),S),']')
  vals <- matrix(vals,nrow=S)
  diag(vals) <- NA
  phi_mult <- c(na.omit(as.vector(vals)))

  #Cond binom
  phi_cb <- c(paste0('phi[',0:(S-1),']'),paste0('R[',0:(S-1),']'))

  umf@phiOrder <- list(multinomial=phi_mult,cond_binom=phi_cb)

  umf
}

#Construct a unmarkedFrameOccuTTD object (time-to-detection occupancy)
unmarkedFrameOccuTTD <- function(y, surveyLength, siteCovs = NULL, obsCovs = NULL,
                                 numPrimary = 1, yearlySiteCovs = NULL)
{

  #Convert vector y to 1 col matrix
  if(is.null(dim(y))) y <- matrix(y, ncol=1)

  #Check if surveyLength input is correct
  if(length(surveyLength)>1 & is.null(dim(surveyLength))){
    surveyLength <- matrix(surveyLength, ncol=1)
  }

  if(length(surveyLength)==1){
    surveyLength <- matrix(surveyLength, nrow=nrow(y), ncol=ncol(y))
  } else if(!identical(dim(surveyLength), dim(y))){
    stop("surveyLength must be a single value or have same dimensions as y",
         call. = FALSE)
  }

  #Check if numPrimary makes sense
  if((ncol(y) %% numPrimary) != 0){
    stop(paste("Uneven number of secondary periods in each primary period.",
               "You can pad with columns of NAs if necessary."),call.=FALSE)
  }

  #Build umf
  umf <- unmarkedMultFrame(y, siteCovs, obsCovs, numPrimary, yearlySiteCovs)
  umf <- as(umf, "unmarkedFrameOccuTTD")
  umf@surveyLength <- surveyLength

  umf
}


# This function constructs an unmarkedMultFrame object.
unmarkedFrameGMM <- function(y, siteCovs = NULL, obsCovs = NULL, numPrimary,
	yearlySiteCovs = NULL, type, obsToY, piFun)
{
    J <- ncol(y) / numPrimary
    if(!missing(type)) {
      if(!type %in% c("removal", "double", "depDouble"))
        stop("if specifying type, it should either be 'removal', 'double', or 'depDouble'")
      switch(type,
        removal = {
          obsToY <- diag(J)
          obsToY[upper.tri(obsToY)] <- 1
          obsToY <- kronecker(diag(numPrimary), obsToY)
          piFun <- "removalPiFun"
          },
        double = {
          obsToY <- matrix(1, 2, 3)
          obsToY <- kronecker(diag(numPrimary), obsToY)
          piFun <- "doublePiFun"
          },
        depDouble = {
          if(J!=2){
            stop("y must have exactly 2 columns per primary period", call.=FALSE)
          }
          obsToY <- matrix(1, 2, 2)
          obsToY <- kronecker(diag(numPrimary), obsToY)
          piFun <- "depDoublePiFun"
          })
    } else {
       ### 2/4/2019 ... cat("this is a test message",fill=TRUE)
       type <- "userDefined"
        # ncol(obsToY) should be ncol(y) (numPrimary*J)
        if(missing(obsToY))
            stop("obsToY is required for gmultmix data with no specified type.")
        if(ncol(obsToY) == J){
           obsToY <- kronecker(diag(numPrimary), obsToY)
       }
        if(ncol(obsToY) == ncol(y)){
           # this is a good situation
        }
        if(ncol(obsToY) != ncol(y) &  ncol(obsToY) != J){
           stop("obsToY wrong dimension.")
       }

      }

    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    umf@yearlySiteCovs <- covsToDF(yearlySiteCovs, "yearlySiteCovs",
                                   numPrimary, nrow(y))
    umf <- as(umf, "unmarkedFrameGMM")
    umf@piFun <- piFun
    umf@samplingMethod <- type
    umf <- umf_to_factor(umf)
    umf
}



# This function constructs an unmarkedMultFrame object.
unmarkedFrameDSO <- function(y, siteCovs=NULL, yearlySiteCovs=NULL, numPrimary,
                             primaryPeriod, dist.breaks, tlength, survey,
                             unitsIn)
{
    J <- ncol(y) / numPrimary
    obsToY <- matrix(1, 1, J)
    obsToY <- kronecker(diag(numPrimary), obsToY)
    if(missing(siteCovs))
        siteCovs <- NULL
    M <- nrow(y)
    T <- numPrimary

    if(missing(primaryPeriod))
        primaryPeriod <- matrix(1:T, M, T, byrow=TRUE)
    if(nrow(primaryPeriod) != M | ncol(primaryPeriod) != T)
        stop("Dimensions of primaryPeriod matrix should be MxT")
    if(any(primaryPeriod < 0, na.rm=TRUE))
        stop("Negative primaryPeriod values are not allowed.")
    if(any(is.na(primaryPeriod)))
        stop("Missing values are not allowed in primaryPeriod.")
    if(!identical(typeof(primaryPeriod), "integer")) {
        mode(primaryPeriod) <- "integer"
        warning("primaryPeriod values have been converted to integers")
        }

    #######well fuck, not sure why this needs to be this way:
    #obsToY = diag(J*T)
    ### This is from unmarkedFramePCO

    umf <- unmarkedFrame(y = y, siteCovs = siteCovs, obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(inherits(yearlySiteCovs, "list")) {
        yearlySiteVars <- names(yearlySiteCovs)
        for(i in seq(length(yearlySiteVars))) {
            if(!inherits(yearlySiteCovs[[i]], c("matrix","data.frame")))
                stop("At least one element of yearlySiteCovs is not a matrix or data frame.")
            if(ncol(yearlySiteCovs[[i]]) != numPrimary |
                nrow(yearlySiteCovs[[i]]) != nrow(y))
                    stop("At least one matrix in yearlySiteCovs has incorrect number of dimensions.")
            }
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs, function(x)
            as.vector(t(x))))
        }
    if(identical(survey, "point")) {
        if(!missing(tlength))
            stop("tlength cannot be specified with point transect data")
        tlength <- rep(1, nrow(y))
        }

    if(identical(survey, "line")) {
      if(length(tlength) != nrow(y)) {
        stop("tlength should be a vector with length(tlength)==nrow(y)")
      }
    }

    if(length(dist.breaks) != J+1){
      stop(paste("dist.breaks should have",J+1,"values"))
    }

    umf@yearlySiteCovs <- yearlySiteCovs
    umf <- as(umf, "unmarkedFrameDSO")
    umf@dist.breaks <- dist.breaks
    umf@survey <- survey
    umf@unitsIn <- unitsIn
    umf@tlength <- tlength
    umf@primaryPeriod <- primaryPeriod
    umf <- umf_to_factor(umf)
    umf
}











unmarkedFrameGDS<- function(y, siteCovs, numPrimary,
	yearlySiteCovs, dist.breaks, survey, unitsIn, tlength)
{
    J <- ncol(y) / numPrimary
    obsToY <- matrix(1, 1, J)
    obsToY <- kronecker(diag(numPrimary), obsToY)
    if(missing(siteCovs))
        siteCovs <- NULL

    umf <- unmarkedFrame(y = y, siteCovs = siteCovs, obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(missing(yearlySiteCovs))
        yearlySiteCovs <- NULL
    umf@yearlySiteCovs <- covsToDF(yearlySiteCovs, "yearlySiteCovs",
                                   numPrimary, nrow(y))

    if(identical(survey, "point")) {
        if(!missing(tlength))
            stop("tlength cannot be specified with point transect data")
        tlength <- rep(1, nrow(y))
        }
    if(identical(survey, "line")) {
      if(length(tlength) != nrow(y)) {
        stop("tlength should be a vector with length(tlength)==nrow(y)")
      }
    }

    if(length(dist.breaks) != J+1){
      stop(paste("dist.breaks should have",J+1,"values"))
    }

    umf <- as(umf, "unmarkedFrameGDS")
    umf@dist.breaks <- dist.breaks
    umf@survey <- survey
    umf@unitsIn <- unitsIn
    umf@tlength <- tlength
    umf <- umf_to_factor(umf)
    umf
}





































# This function constructs an  object.
unmarkedFrameGPC <- function(y, siteCovs=NULL, obsCovs=NULL, numPrimary,
                             yearlySiteCovs=NULL) {
    if(numPrimary < 2)
        stop("numPrimary must be >1. Use pcount of numPrimary=1")
    J <- ncol(y) / numPrimary
    obsToY <- diag(J*numPrimary)
    if(missing(siteCovs))
        siteCovs <- NULL

    umf <- unmarkedFrame(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
                         obsToY = obsToY)
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    if(missing(yearlySiteCovs))
        yearlySiteCovs <- NULL
    umf@yearlySiteCovs <- covsToDF(yearlySiteCovs, "yearlySiteCovs",
                                   numPrimary, nrow(y))
    umf <- as(umf, "unmarkedFrameGPC")
    umf <- umf_to_factor(umf)
    umf
}





unmarkedFramePCO <- function(y, siteCovs = NULL, obsCovs = NULL,
    yearlySiteCovs = NULL, mapInfo, numPrimary, primaryPeriod)
{
    M <- nrow(y)
    T <- numPrimary
    J <- ncol(y) / T
    if(missing(primaryPeriod))
        primaryPeriod <- matrix(1:T, M, T, byrow=TRUE)
    if(nrow(primaryPeriod) != M | ncol(primaryPeriod) != T)
        stop("Dimensions of primaryPeriod matrix should be MxT")
    if(any(primaryPeriod < 0, na.rm=TRUE))
        stop("Negative primaryPeriod values are not allowed.")
    if(any(is.na(primaryPeriod)))
        stop("Missing values are not allowed in primaryPeriod.")
    if(!identical(typeof(primaryPeriod), "integer")) {
        mode(primaryPeriod) <- "integer"
        warning("primaryPeriod values have been converted to integers")
        }
    ya <- array(y, c(M, J, T))
    yt.na <- apply(!is.na(ya), c(1,3), any)
    yt.na <- which(!yt.na)
    d.na <- which(is.na(primaryPeriod))
    if(!all(d.na %in% yt.na))
        stop("primaryPeriod values must be supplied for all non-missing values of y")
    increasing <- function(x) {
        x <- x[!is.na(x)]
        all(order(x) == 1:length(x))
        }
    if(!all(apply(primaryPeriod, 1, increasing)))
        stop("primaryPeriod values must increase over time for each site")

    obsCovs <- covsToDF(obsCovs, "obsCovs", J*T, M)
    umf <- unmarkedFrame(y, siteCovs, obsCovs, obsToY = diag(J*T))
    umf <- as(umf, "unmarkedMultFrame")
    umf@numPrimary <- numPrimary
    umf@yearlySiteCovs <- covsToDF(yearlySiteCovs, "yearlySiteCovs",
                                   numPrimary, nrow(y))
    umf <- as(umf, "unmarkedFramePCO")
    umf@primaryPeriod <- primaryPeriod
    # There is no obsToY function
    umf <- umf_to_factor(umf)
    return(umf)
}


unmarkedFrameMMO <- function(y, siteCovs = NULL, obsCovs = NULL,
    yearlySiteCovs = NULL, numPrimary, type, primaryPeriod)
{

    M <- nrow(y)
    T <- numPrimary
    J <- ncol(y) / T

    if(missing(primaryPeriod))
        primaryPeriod <- matrix(1:T, M, T, byrow=TRUE)
    if(nrow(primaryPeriod) != M | ncol(primaryPeriod) != T)
        stop("Dimensions of primaryPeriod matrix should be MxT")
    if(any(primaryPeriod < 0, na.rm=TRUE))
        stop("Negative primaryPeriod values are not allowed.")
    if(any(is.na(primaryPeriod)))
        stop("Missing values are not allowed in primaryPeriod.")
    if(!identical(typeof(primaryPeriod), "integer")) {
        mode(primaryPeriod) <- "integer"
        warning("primaryPeriod values have been converted to integers")
        }

    ya <- array(y, c(M, J, T))
    yt.na <- apply(!is.na(ya), c(1,3), any)
    yt.na <- which(!yt.na)
    d.na <- which(is.na(primaryPeriod))
    if(!all(d.na %in% yt.na))
        stop("primaryPeriod values must be supplied for all non-missing values of y")
    increasing <- function(x) {
        x <- x[!is.na(x)]
        all(order(x) == 1:length(x))
        }
    if(!all(apply(primaryPeriod, 1, increasing)))
        stop("primaryPeriod values must increase over time for each site")

    umf <- unmarkedFrameGMM(y, siteCovs, obsCovs, numPrimary, yearlySiteCovs,
                            type)
    umf <- as(umf, "unmarkedFrameMMO")
    umf@primaryPeriod <- primaryPeriod

    umf
}



################ SHOW METHODS ############################################


setMethod("show", "unmarkedFrame", function(object)
{
    df <- as(object, "data.frame")
    cat("Data frame representation of unmarkedFrame object.\n")
    print(df)
})


setMethod("show", "unmarkedMultFrame",
    function(object)
{
    df <- as(object, "data.frame")
    cat("Data frame representation of unmarkedFrame object.\n")
    print(df)
})

setMethod("show", "unmarkedFrameOccuMulti", function(object)
{
    df <- as(object, "data.frame")
    cat("Data frame representation of unmarkedFrame object.\n")
    cat("Only showing observation matrix for species 1.\n")
    print(df)
})

setMethod("show", "unmarkedFrameOccuTTD", function(object)
{
  df <- as(object, "data.frame")
  sl <- object@surveyLength
  if(ncol(sl)>1){
    colnames(sl) <- paste0('surveyLength',1:ncol(sl))
  } else {
    colnames(sl) <- 'surveyLength'
  }
  df <- cbind(df, sl)
  cat("Data frame representation of unmarkedFrame object.\n")
  print(df)
})

############################ EXTRACTORS ##################################

# Extractor for site level covariates
setGeneric("siteCovs", function(object,...) standardGeneric("siteCovs"))

setMethod("siteCovs", "unmarkedFrame", function(object) {
    return(object@siteCovs)
})

setGeneric("yearlySiteCovs", function(object,...)
    standardGeneric("yearlySiteCovs"))
setMethod("yearlySiteCovs", "unmarkedMultFrame", function(object) {
    return(object@yearlySiteCovs)
})


setGeneric("obsCovs", function(object,...) standardGeneric("obsCovs"))
setMethod("obsCovs", "unmarkedFrame", function(object, matrices = FALSE) {
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


setGeneric("obsCovs<-", function(object, value)
    standardGeneric("obsCovs<-"))
setReplaceMethod("obsCovs", "unmarkedFrame", function(object, value) {
    if(inherits(object, "unmarkedFrameDS"))
        stop("unmarkedFrameDS objects cannot have obsCovs")
    object@obsCovs <- as.data.frame(value)
    object
})


setGeneric("siteCovs<-", function(object, value)
    standardGeneric("siteCovs<-"))
setReplaceMethod("siteCovs", "unmarkedFrame", function(object, value) {
    object@siteCovs <- as.data.frame(value)
    object
})


setGeneric("yearlySiteCovs<-",
	function(object, value) standardGeneric("yearlySiteCovs<-"))
setReplaceMethod("yearlySiteCovs", "unmarkedMultFrame",
    function(object, value) {
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
setMethod("coordinates", "unmarkedFrame", function(object) {
    object@mapInfo@coordinates
})


setGeneric("projection", function(object) standardGeneric("projection"))
setMethod("projection", "unmarkedFrame", function(object) {
    object@mapInfo@projection
})

################################### SUMMARY METHODS ######################


setMethod("summary", "unmarkedFrame", function(object,...) {
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



setMethod("summary", "unmarkedFrameDS", function(object, ...)
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




setMethod("summary", "unmarkedMultFrame", function(object,...) {
    cat("unmarkedFrame Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat("Maximum number of observations per site:",ncol(object@y),"\n")
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


setMethod("summary", "unmarkedFrameOccuMulti", function(object,...) {
    cat("unmarkedFrame Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat(length(object@ylist), "species:", names(object@ylist),"\n")
    cat("Maximum number of observations per site:",obsNum(object),"\n")
    mean.obs <- sapply(object@ylist,function(x) mean(rowSums(!is.na(x))))
    cat("Mean number of observations per site:\n")
    invisible(sapply(1:length(mean.obs), function(x)
        cat(paste(names(mean.obs)[x],': ',mean.obs[x],'  ',sep='')))); cat("\n")
    s.one <- sapply(object@ylist,function(x) sum(rowSums(x,na.rm=T)>0))
    cat("Sites with at least one detection:\n")
    invisible(sapply(1:length(s.one), function(x)
        cat(paste(names(s.one)[x],': ',s.one[x],'  ',sep='')))); cat("\n")
    cat("Tabulation of y observations:\n")
    for (i in 1:length(object@ylist)){
      cat(paste(names(object@ylist)[i],':',sep=''))
      print(table(object@ylist[[i]], exclude=NULL))
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

setMethod("summary", "unmarkedFrameOccuTTD", function(object,...) {

    cat("unmarkedFrame Object\n\n")
    cat(nrow(object@y), "sites\n")
    cat("Maximum number of surveys per site:",ncol(object@y),"\n")
    mean.obs <- mean(rowSums(!is.na(getY(object))))
    cat("Mean number of surveys per site:",round(mean.obs,2),"\n")
    cat("Number of primary survey periods:", object@numPrimary, "\n")
    cat("Number of secondary surveys:",
        obsNum(object) / object@numPrimary, "\n")

    detect <- getY(object)<object@surveyLength
    cat("Sites with at least one detection:",
        sum(apply(detect, 1, function(x) any(x, na.rm=TRUE))),
        "\n\n")
    cat("Distribution of detection times:\n")
    detect_times <- as.numeric(getY(object)[detect])
    print(summary(detect_times, exclude=NULL))
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


################################# PLOT METHODS ###########################
# TODO:  come up with nice show/summary/plot methods for each data types.

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
    y2 <- reshape(y, idvar=c("site", "group"), varying=list(1:ncol(getY(x))),
              v.names="value", direction="long")
    y2$variable <- factor(paste("obs", y2$time), levels=colnames(y))
    if(missing(colorkey))
        colorkey <- list(at=0:(ym+1), labels=list(labels=as.character(0:ym),
            at=(0:ym)+0.5))
    levelplot(value ~ variable*site | group, y2,
        scales=list(relation="free", x=list(labels=1:J)),
        colorkey=colorkey, strip=strip, xlab=xlab, ylab=ylab, ...)
})

setMethod("plot", c(x="unmarkedFrameOccuMulti", y="missing"),
	function (x, y, colorkey, ylab="Site", xlab="Observation", ...)
{
    y <- getY(x)
    ym <- max(y, na.rm=TRUE)
    M <- nrow(y)
    J <- ncol(y)
    S <- length(x@ylist)
    y <- as.data.frame(do.call(rbind,x@ylist))
    colnames(y) <- paste("obs",1:J)
    y$site <- rep(1:M,S)
    y$species <- as.factor(rep(names(x@ylist),each=M))
    y2 <- reshape(y, idvar=c("site", "species"), varying=list(1:obsNum(x)),
              v.names="value", direction="long")
    y2$variable <- factor(paste("obs", y2$time), levels=colnames(y))
    if(missing(colorkey))
        colorkey <- list(at=0:(ym+1), labels=list(labels=as.character(0:ym),
            at=(0:ym)+0.5))
    levelplot(value ~ variable*site | species, y2,
        scales=list(relation="free", x=list(labels=1:J)),
        colorkey=colorkey, strip=T, xlab=xlab, ylab=ylab,
        labels=names(x@ylist), ...)
})

setMethod("plot", c("unmarkedFrameOccuTTD", y="missing"),
  function(x, y, ...){

  y <- getY(x)
  y <- y[y<x@surveyLength]
  hist(y, xlab="Time to first detection", ylab="Frequency", main="", ...)
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



################################# SELECTORS ##############################

# i is the vector of sites to extract

setMethod("[", c("unmarkedFrame", "numeric", "missing", "missing"),
    function(x, i)
{
    M <- numSites(x)
    if(length(i) == 0) return(x)
    if(any(i < 0) && any(i > 0))
        stop("i must be all positive or all negative indices.")
    if(all(i < 0)) { # if i is negative, then convert to positive
        i <- (1:M)[i]
        }
    y <- getY(x)[i,,drop=FALSE]
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


## remove obs only
### RBC: Why??? this doesn't allow umf[,c(1,1)]
setMethod("[", c("unmarkedFrame", "missing", "numeric", "missing"),
		function(x, i, j)
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


# i is as before and j is the obsNum to remove and corresponding y's
setMethod("[", c("unmarkedFrame","numeric", "numeric", "missing"),
		function(x, i, j)
{
    ## first remove sites
    umf <- x[i,]
    umf <- umf[,j]
    umf
})



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


#[ Methods for multispecies occupancy frames
setMethod("[", c("unmarkedFrameOccuMulti", "numeric", "missing", "missing"),
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

setMethod("[", c("unmarkedFrameOccuMulti", "missing", "numeric", "missing"),
		function(x, i, j)
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
setMethod("[", c("unmarkedMultFrame", "missing", "numeric", "missing"),
		function(x, i, j)
{
    J <- obsNum(x)/x@numPrimary
    obs <- rep(1:x@numPrimary, each = J)
    years <- 1:x@numPrimary
    numPrimary <- length(j)
    obsj <- match(obs, j)
    j2 <- which(!is.na(obsj))
    u <- callNextMethod(x, i, j2)
    ysc <- yearlySiteCovs(x)
    if(!is.null(ysc)) {
        ysc <- ysc[rep(!is.na(match(years, j)), nrow(getY(x))),, drop=FALSE]
        u@yearlySiteCovs <- ysc
        }
    u@numPrimary <- numPrimary
    return(u)
})



## for multframes, must remove years at a time
setMethod("[", c("unmarkedMultFrame", "numeric", "missing", "missing"),
		function(x, i, j)
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


setMethod("[", c("unmarkedFrameOccuMS", "numeric", "missing", "missing"),
    function(x, i, j)
{
  multf <- callNextMethod(x, i, j)
  unmarkedFrameOccuMS(y=getY(multf), siteCovs=siteCovs(multf),
                      yearlySiteCovs=yearlySiteCovs(multf),
                      obsCovs=obsCovs(multf),
                      numPrimary=x@numPrimary)
})

setMethod("[", c("unmarkedFrameGMM", "numeric", "missing", "missing"),
		function(x, i, j)
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


setMethod("[", c("unmarkedFrameGPC", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    class(multf) <- "unmarkedFrameGPC"
    multf
})


setMethod("[", c("unmarkedFrameGPC", "missing", "numeric", "missing"),
		function(x, i, j)
{
    multf <- as(x, "unmarkedMultFrame")
    out <- callNextMethod(multf, i, j) # unmarkedMultFrame
    as(out, "unmarkedFrameGPC")
})





setMethod("[", c("unmarkedFrameGDS", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
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



setMethod("[", c("unmarkedFramePCO", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    unmarkedFramePCO(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=x@numPrimary,
                     primaryPeriod=x@primaryPeriod[i,,drop=FALSE])
})


setMethod("[", c("unmarkedFramePCO", "missing", "numeric", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    unmarkedFramePCO(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=length(j),
                     primaryPeriod=x@primaryPeriod[,j,drop=FALSE])
})


setMethod("[", c("unmarkedFrameOccuTTD", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
    unmarkedFrameOccuTTD(y=getY(multf), siteCovs=siteCovs(multf),
                     yearlySiteCovs=yearlySiteCovs(multf),
                     obsCovs=obsCovs(multf),
                     numPrimary=x@numPrimary,
                     surveyLength=x@surveyLength[i,,drop=FALSE])
})


setMethod("[", c("unmarkedFrameDSO", "numeric", "missing", "missing"),
		function(x, i, j)
{
    multf <- callNextMethod(x, i, j) # unmarkedMultFrame
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


setMethod("[", c("unmarkedFrameOccuTTD", "missing", "numeric", "missing"),
		function(x, i, j)
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




setMethod("head", "unmarkedFrame", function(x, n) {
    if(missing(n)) n <- 10
    umf <- x[1:n,]
    umf
})

############################### COERCION #################################

setAs("data.frame", "unmarkedFrame", function(from)
{
    umf <- formatWide(from)
    umf
})



setAs("unmarkedFrame", "data.frame", function(from)
{
    obsCovs <- obsCovs(from)
    siteCovs <- siteCovs(from)
    if(inherits(from, 'unmarkedMultFrame')) {
      yearlySiteCovs <- yearlySiteCovs(from)
    } else {
      yearlySiteCovs <- NULL
    }
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
    if(is.null(yearlySiteCovs)) {
        yearlySiteCovs <- matrix(0,nrow(y),0)
    } else {
        yearlySiteCovs <- data.frame(lapply(yearlySiteCovs,
            function(x) matrix(x, nrow(y), length(x)/nrow(y), byrow=T)))
        }
    df <- data.frame(y, siteCovs, yearlySiteCovs, obsCovs)
    df
})






