
genFixedNLL <- function(nll, whichFixed, fixedValues)
{
    function(params) {
        params[whichFixed] <- fixedValues
        do.call(nll, list(params))
        }
}

# nll the original negative log likelihood function
# MLE the full vector of MLE values
profileCI <- function(nll, whichPar, MLE, interval, level)
{
    stopifnot(length(whichPar) == 1)
    MLEnll <- nll(MLE)
    nPar <- length(MLE)
	chsq <- qchisq(level, 1)/2
    f <- function(value) {
        fixedNLL <- genFixedNLL(nll, whichPar, value)
            mleRestricted <- optim(MLE, fixedNLL)$value
        mleRestricted - MLEnll - chsq
        }
    lower <- tryCatch(uniroot(f, c(interval[1],MLE[whichPar]))$root,
        error = function(e) {
            warning("Lower endpoint of profile confidence interval is on the boundary.",
        call. = FALSE)
        -Inf
        })
    upper <- tryCatch(upper <- uniroot(f, c(MLE[whichPar], interval[2]))$root,
        error = function(e) {
            warning("Upper endpoint of profile confidence interval is on the boundary.",
        call. = FALSE)
        Inf
        })

    return(c(lower,upper))
}

## link functions and their gradients
logistic <- function(x) {
  1/(1 + exp(-x))
}


logistic.grad <- function(x) {
  exp(-x)/(exp(-x)+1)^2
}


log.grad <- function(x) { # duh! (but for clarity)
  1/x
}


explink <- function(x) exp(x)


identLink <- function(x) x


identLinkGrad <- function(x) 1

# !!! This is now done in C++. See rowProds() below
## use logarithms to vectorize row-wise products
## this speeds things up a LOT (vs. apply(x,1,prod))

#rowProds <-
#function(x, na.rm = FALSE)
#{
#  exp(rowSums(log(x), na.rm = na.rm))
#}

# helper function to coerce an array of matrices to a list

arrToList <- function(x){
  nl <- list()
  for(i in 1:dim(x)[3]) {
    nl[[i]] <- x[,,i]
  }
  names(nl) <- dimnames(x)[[3]]
  nl
}

## compute estimated asymptotic variances of parameter estimates
## using the observed information matrix

#sd.est <- function(fm) {
#    sqrt(diag(solve(fm$hessian)))
#}

## delta method for variance of proportion given variance of its logistic-
## transformed counterpart
##' @nord
#sd.prop <- function(est,sd.est) {
#    exp(-est)/(1 + exp(-est))^2 * sd.est
#}

### track linked list of parameters using a data frame
### add row to linked list

addParm <- function(list.df, parm.name, parm.length) {
    if(parm.length > 0) {
        if(nrow(list.df) == 0) {
            last.ind <- 0
        } else {
            last.ind <- list.df$end[nrow(list.df)]
        }
        parm.df <- data.frame(parameter = parm.name, start = last.ind + 1,
                              end = last.ind + parm.length,
                              stringsAsFactors = FALSE)
        list.df <- rbind(list.df, parm.df)
    }
    return(list.df)
}


parmNames <- function(list.df) {
    npar <- list.df$end[nrow(list.df)]
    names <- character(npar)
    for(i in 1:npar) {
        which.par <- which(i >= list.df$start & i <= list.df$end)
        names[i] <- list.df$parameter[which.par]
    }
    return(names)
}


# This function converts an appropriatedly formated comma-separated
# values file (.csv) to a format usable by \emph{unmarked}'s fitting
# functions (see \emph{Details}).
csvToUMF <-
function(filename, long=FALSE, type, species = NULL, ...)
{
  dfin <- read.csv(filename)

  if(long == TRUE) return(formatLong(dfin, species, type = type, ...))
  else return(formatWide(dfin, type = type, ...))
}

# utility function to create a variable that follows the dates as 1,2,3,...
# site id is first column
# julian date is second column

dateToObs <-
function(dfin)
{
  sitecol <- dfin[[1]]
  datecol <- dfin[[2]]


  # order by site, then obs date
  dfin <- dfin[order(sitecol,datecol),]
  sitecol <- dfin[[1]]
  datecol <- dfin[[2]]

  dTab <- table(datecol,sitecol)
  sites <- unique(sitecol)
  nSite <- length(sites)
  nStop <- colSums(dTab)
  nStop <- nStop[nStop > 0]  # get rid of the stops for sites with no stops

  obsNum <- numeric(length(sitecol))
  # for each site i, replace stops with 1:nStop[i]
  for(i in 1:nSite){
    stops <- which(sitecol == sites[i])
    obsNum[stops] <- 1:nStop[i]
  }

  dfout <- cbind(dfin,obsNum)
  dfout
}

# take long data set and return data list
# column names must be
# site names, first
# date, one column
# response, one column
# obs vars, one per column
formatLong <-
function(dfin, species = NULL, type)
{

  if(missing(type)) stop("type must be supplied")

  ## copy dates to last column so that they are also a covdata var
  nc <- ncol(dfin)
  dfin[[nc+1]] <- dfin[[2]]
  names(dfin)[nc+1] <- "Date"

  if(!is.null(species)) {
    dfin$y <- ifelse(dfin$species == species, dfin$y, 0)
    dfin$y[is.na(dfin$y)] <- 0
    dfin$species = NULL
  }
# TODO: double check that multiple cells per site*time are handled correctly.
#  # sum up counts within time/site
#  expr <- substitute(recast(dfin[,1:3], sv + dv ~ ..., id.var = 1:2,
#                            fun.aggregate = sum),
#                     list(sv = as.name(names(dfin)[1]),
#                          dv = as.name(names(dfin)[2])))
#  dfin2 <- eval(expr)
#  dfin1 <- dfin[!duplicated(dfin[,1:2]),]
#
#  dfin <- merge(dfin1,dfin2, by = 1:2)
#  dfin[,3] <- dfin[,length(dfin)]
#  dfin <- dfin[,-length(dfin)]
  names(dfin)[3] <- "y"

  dfin <- dateToObs(dfin)
  dfnm <- colnames(dfin)
  nV <- length(dfnm) - 1  # last variable is obsNum
  expr <- substitute(recast(dfin, newvar ~ obsNum + variable,
    id.var = c(dfnm[1],"obsNum"), measure.var = dfnm[3]),
                                list(newvar=as.name(dfnm[1])))
  y <- as.matrix(eval(expr)[,-1])
  attr(y,"class") <- "matrix"

  expr <- substitute(recast(dfin, newvar ~ obsNum ~ variable,
    id.var = c(dfnm[1],"obsNum"), measure.var = dfnm[4:nV]),
                                list(newvar=as.name(dfnm[1])))
  obsvars <- eval(expr)
  which.date <- which(dimnames(obsvars)$variable == "Date")
  dimnames(obsvars)$variable[which.date] <- "JulianDate"

  obsvars.matlist <- arrToList(obsvars)
  obsvars.veclist <- lapply(obsvars.matlist, function(x) as.vector(t(x)))
  obsvars.df <- data.frame(obsvars.veclist)

	do.call(type, list(y = y, obsCovs = obsvars.df))
}

# column names must be
# site (optional, but if present, labeled "site")
# response: y.1, y.2, ..., y.J
# site vars: namefoo, namebar, ...
# obs vars: namefoo.1, namefoo.2, ..., namefoo.J, namebar.1, ..., namebar.J,...

formatWide <-
function(dfin, sep = ".", obsToY, type, ...)
{
	# escape separater if it is regexp special
	reg.specials <- c('.', '\\', ':', '|', '(', ')', '[', '{', '^', '$', '*', '+', '?')
	if(sep %in% reg.specials) {
		sep.reg <- paste("\\",sep,sep="")
	} else {
		sep.reg <- sep
	}

	dfnm <- colnames(dfin)

	y <- grep(paste("^y",sep.reg,"[[:digit:]]", sep=""),dfnm)
	J <- length(y)
	y <- as.matrix(dfin[,y])
	M <- nrow(y)

  if(identical(tolower(colnames(dfin))[1],"site")) {
		dfin <- dfin[,-1]
		dfnm <- dfnm[-1]
	}

  ncols <- length(dfnm)
	obsCovsPresent <- FALSE
	siteCovsPresent <- FALSE
  i <- J + 1
  while(i <= ncols) {     # loop through columns
    if(!identical(grep(paste(sep.reg,"[[:digit:]]+$",sep=""),dfnm[i]),integer(0))) {  # check if this is obsdata
      newvar.name <- sub(paste(sep.reg,"[[:digit:]]+$",sep=""),'',dfnm[i])
			newvar <- dfin[,grep(paste(newvar.name,sep.reg,"[[:digit:]]+$",sep=""),dfnm)]
			if(obsCovsPresent) {
					if(ncol(newvar) != R) {
						stop("Not all observation-level covariates have the same number of columns.")
					} else {
						obsCovs[newvar.name] <- as.vector(t(newvar))
					}
			} else {
				obsCovsPresent <- TRUE
				R <- ncol(newvar)
				obsCovs <- data.frame(newvar = as.vector(t(newvar)))
			}
			colnames(obsCovs)[length(obsCovs)] <- newvar.name
      i <- i + R
    }
    else {
			if(siteCovsPresent){
				siteCovs <- cbind(siteCovs,dfin[,i])
			} else {
				siteCovsPresent <- TRUE
				siteCovs <- data.frame(newvar = dfin[,i])
			}
      colnames(siteCovs)[length(siteCovs)] <- dfnm[i]
      i <- i + 1
    }
  }

	if(!obsCovsPresent) obsCovs <- NULL
	if(!siteCovsPresent) siteCovs <- NULL

	## if we don't know the true obsToY yet, use R X J matrix of ones or diag if R == J
	if(missing(obsToY)) {
		if(identical(R,J)) {
			obsToY <- diag(J)
		} else {
			obsToY <- matrix(0, R, J)
		}
	}

	do.call(type, list(y = y, siteCovs = siteCovs, obsCovs = obsCovs, ...))
}


# This convenience function converts multi-year data in long format to unmarkedMultFrame Object.  See Details for more information.

formatMult <-
function(df.in)
{
  years <- sort(unique(df.in[[1]]))
  nY <- length(years)
  df.obs <- list()
  nsamp <- numeric()
  maxsamp <- max(table(df.in[[1]], df.in[[2]])) # the maximum samples/yr
  for(t in 1:nY){
    df.t <- df.in[df.in[[1]] == years[t],] # subset for current year
    df.t <- df.t[,-1] # remove year column
    df.t <- dateToObs(df.t)
    nsamp <- max(df.t$obsNum)
    if(nsamp < maxsamp) {
      newrows <- df.t[1:(maxsamp - nsamp), ] # just a placeholder
      newrows[,"obsNum"] <- ((nsamp + 1) : maxsamp)
      newrows[,3 : (ncol(df.t) - 1)] <- NA
      df.t <- rbind(df.t, newrows)
    }
    df.obs <- rbind(df.obs,cbind(year = years[t],df.t))
  }
  dfnm <- colnames(df.obs)
  nV <- length(dfnm) - 1  # last variable is obsNum

  # create y matrix using reshape
  expr <- substitute(recast(df.obs, var1 ~ year + obsNum + variable,
    id.var = c(dfnm[2],"year","obsNum"), measure.var = dfnm[4]),
                                list(var1 = as.name(dfnm[2])))
  y <- as.matrix(eval(expr)[,-1])

  # create obsdata with reshape
  # include date (3rd col) and other measured vars
  expr <- substitute(recast(df.obs, newvar ~ year + obsNum ~ variable,
    id.var = c(dfnm[2],"year","obsNum"), measure.var = dfnm[c(3,5:nV)]),
                                list(newvar=as.name(dfnm[2])))
  obsvars <- eval(expr)

  rownames(y) <- dimnames(obsvars)[[1]]
  colnames(y) <- dimnames(obsvars)[[2]]
	y <- as.matrix(y)

	obsvars.list <- arrToList(obsvars)
	obsvars.list <- lapply(obsvars.list, function(x) as.vector(t(x)))
	obsvars.df <- as.data.frame(obsvars.list)

	## check for siteCovs
	obsNum <- ncol(y)
	M <- nrow(y)
	site.inds <- matrix(1:(M*obsNum), M, obsNum, byrow = TRUE)
	siteCovs <- sapply(obsvars.df, function(x) {
				obsmat <- matrix(x, M, obsNum, byrow = TRUE)
				l.u <- apply(obsmat, 1, function(y) {
							row.u <- unique(y)
							length(row.u[!is.na(row.u)])
						})
				if(all(l.u %in% 0:1)) {  ## if there are 0 or 1 unique vals per row, we have a sitecov
					u <- apply(obsmat, 1, function(y) {
								row.u <- unique(y)
								if(!all(is.na(row.u)))  ## only remove NAs if there are some non-NAs.
									row.u <- row.u[!is.na(row.u)]
								row.u
							})
					u
				}
			})
	siteCovs <- as.data.frame(siteCovs[!sapply(siteCovs, is.null)])
	if(nrow(siteCovs) == 0) siteCovs <- NULL

  ## only check non-sitecovs
  obsvars.df2 <- as.data.frame(obsvars.df[, !(names(obsvars.df) %in% names(siteCovs))])
  names(obsvars.df2) <- names(obsvars.df)[!(names(obsvars.df) %in% names(siteCovs))]

	yearlySiteCovs <- sapply(obsvars.df2, function(x) {
				obsmat <- matrix(x, M*nY, obsNum/nY, byrow = TRUE)
				l.u <- apply(obsmat, 1, function(y) {
							row.u <- unique(y)
							length(row.u[!is.na(row.u)])
						})
				if(all(l.u %in% 0:1)) {  ## if there are 0 or 1 unique vals per row, we have a sitecov
					u <- apply(obsmat, 1, function(y) {
								row.u <- unique(y)
								if(!all(is.na(row.u)))  ## only remove NAs if there are some non-NAs.
									row.u <- row.u[!is.na(row.u)]
								row.u
							})
					u
				}
			})
	yearlySiteCovs <- as.data.frame(yearlySiteCovs[!sapply(yearlySiteCovs, is.null)])
	if(nrow(yearlySiteCovs) == 0) yearlySiteCovs <- NULL

	umf <- unmarkedMultFrame(y = y, siteCovs = siteCovs, obsCovs = obsvars.df, yearlySiteCovs = yearlySiteCovs,
			numPrimary = nY)
  return(umf)
}

# function to take data of form
# site  | species | count
# to
# site | spp1 | spp2 | ...

sppLongToWide <-
function(df.in)
{
  df.m <- melt(df.in, id = c("site", "spp"))
  df.out <- cast(df.m, site ~ spp, add.missing=T, fill = 0)
  df.out <- df.out[order(df.out$site),]
  df.out
}

# get estimated psi from rn fit

getPsi <-
function(lam)
{
  1-exp(-lam)
}

# get estimatd p from rn fit (only for a null type model so far)

getP.bar <-
function(lam, r)
{
  K = 30
  psi <- getPsi(lam)
  pN.k <- dpois(0:K,lam)
  pY.k <- 1 - (1 - r)^(0:30)
  sum(pY.k * pN.k)
}




## function to move data around:
## converts array obsdata to a list
## copies site covariate info from obsdata to sitedata
## puts all site covariates back into obsdata
## needed because all site vars need to be available for both state and det models
##' @nord
#arrangeData <-
#function(data)
#{
#  require(abind)
#  y <- data$y
#  sitedata <- data$covdata.site
#  obsdata <- data$covdata.obs
#
#  J <- ncol(y)
#  M <- nrow(y)
#  nSV <- length(sitedata)
#
#  # if not null, then just add "ones"
#  if(!is.null(obsdata)) {
#      if(class(obsdata) == "list") obsdata$ones <- matrix(1,M,J)
#      if(class(obsdata) == "array") {
#          obsdata <- abind(obsdata, ones = matrix(1,M,J))
#      }
#  }
#  if(!is.null(sitedata)) sitedata <- cbind(ones = rep(1,M),sitedata)
#
#  # if data components are null, create as just ones
#  if(is.null(obsdata)) obsdata <- list(ones = matrix(1,M,J))
#  if(is.null(sitedata)) sitedata=data.frame(ones = rep(1,M))
#
#  # if obsdata is an array, coerce it to a list
#  if(identical(class(obsdata),"array")) obsdata <- arrToList(obsdata)
#  nOV <- length(obsdata)
#
#  # move all site data (single vectors and matrices of J repeated vectors)
#  # in obsdata into sitedata
#  toDel <- numeric(0)
#  nuniq <- function(x) length(as.numeric(na.omit(unique(x)))) # lil' helper fun
#  for(i in 1:nOV){
#    # test for equality across rows (matrix site data)
#    eqRow <- as.numeric(apply(as.matrix(obsdata[[i]]), 1, nuniq) == 1)
#    isRep <- as.logical(prod(eqRow)) # make sure all rows are
#    # move into site data if (vector) or (repeated vector as matrix)
#    if((dim(as.matrix(obsdata[[i]]))[2] == 1) || isRep){
#      obsdata[[i]] <- matrix(obsdata[[i]],nrow = M, ncol = J)
#      sitedata <- cbind(sitedata,obsdata[[i]][,1])
#      colnames(sitedata)[length(sitedata)] <- names(obsdata[i])
#      toDel <- c(toDel,i)
#    }
#    # ensure that obsdata is made of matrices rather than dataframes
#    obsdata[[i]] <- as.matrix(obsdata[[i]])
#  }
#  if(length(toDel) > 0) {   #obsdata[[toDel]] <- NULL # remove sitedata from obsdata
#    for(t in toDel) {
#      obsdata[[t]] <- NULL
#    }
#  }
#  if(length(obsdata) == 0) obsdata <- list(ones = matrix(1,M,J))
#  nSV <- length(sitedata) # update nSV
#  nOV <- length(obsdata)
#
#  # make all site terms into obs terms by copying them to
#  # obsdata (from vector to a matrix of repeated vectors)
#  # needed if site variables are used to model detection.
#  for(i in 1:nSV){
#    obsdata[[nOV + i]] <- matrix(sitedata[[i]],nrow=M,ncol=J)
#    names(obsdata)[nOV + i] <- colnames(sitedata)[i]
#  }
#  obsvars <- names(obsdata)
#  nOV <- length(obsdata) # update length
#
#  list(y = y, covdata.site = sitedata, covdata.obs = obsdata)
#}

# Design matrices




meanstate <- function(x) {
    K <- length(x) - 1
    sum(x*(0:K))
}

truncateToBinary <- function(y) {
  if(max(y, na.rm = TRUE) > 1) {
    y <- ifelse(y > 0, 1, 0)
    warning("Some observations were > 1.  These were truncated to 1.")
  }
  return(y)
}


getSS <- function(phi) {
	ev.length <- nrow(phi)
	ev <- tryCatch(eigen(t(phi))$vectors[,1],
			error = function(x) rep(NA, ev.length))
	ev/sum(ev)
}

imputeMissing <- function(umf, whichCovs = seq(length=ncol(obsCovs(umf)))) {


  ## impute observation covariates
  if(!is.null(umf@obsCovs)) {
    obsCovs <- umf@obsCovs
    J <- obsNum(umf)
    M <- nrow(obsCovs)/J
    obs <- as.matrix(obsCovs[,whichCovs])
    whichrows <- apply(obs, 1, function(x) any(!is.na(x)))
    if(sum(whichrows) == 0) return(obsCovs)
    whichels <- matrix(whichrows, M, J, byrow = TRUE)
    for(i in seq(length=length(whichCovs))) {
      obs.i <- obs[,i]
      obs.i.mat <- matrix(obs.i, M, J, byrow = TRUE) # get ith obsvar
      obs.i.missing <- is.na(obs.i.mat) & !whichels
      obs.i.imputed <- obs.i.mat
      for(j in 1:M) {
        for(k in 1:J) {
          if(obs.i.missing[j,k])
            if(all(is.na(obs.i.mat[j,]))) {
              obs.i.imputed[j,k] <- mean(obs.i.mat[,k], na.rm = T)
            } else {
              obs.i.imputed[j,k] <- mean(c(mean(obs.i.mat[j,],na.rm = T),
                                           mean(obs.i.mat[,k], na.rm = T)))
            }
        }
      }
      obsCovs[,whichCovs[i]] <- as.numeric(t(obs.i.imputed))
    }
    umf@obsCovs <- obsCovs
  }
                                        # TODO: impute site covariates
  return(umf)
}


#wideToUMF <- function(formula, data)
#{
#	detformula <- as.formula(formula[[2]])
#	stateformula <- as.formula(paste("~",formula[3],sep=""))
#	yVar <- all.vars(detformula)[1]
#	detVars <- all.vars(detformula)[-1]
#
#	# escape separater if it is regexp special
#	reg.specials <- c('.', '\\', ':', '|', '(', ')', '[', '{', '^', '$', '*', '+', '?')
#	if(sep %in% reg.specials) {
#		sep.reg <- paste("\\",sep,sep="")
#	} else {
#		sep.reg <- sep
#	}
#
#	## y columns must be in format y.1, y.2, ..., y.J where sep="." here
#	yNames <- grep(paste(yVar,sep.reg,"[[:digit:]]+$",sep=""), colnames(data), value=TRUE)
#	timevary <- yNames
#
#	## detection covs must be in format either name alone, or name.1, name.2, ..., name.Q where sep="." here
#	detCols <- lapply(detVars,
#			function(x) sort(grep(paste("^",x,"(",sep.reg,"[[:digit:]]+)?$",sep=""), colnames(data),
#								value=TRUE)))
#
#	names(detCols) <- detVars
#	if (length(detVars) > 0) {
#		toadd <- detVars[sapply(detCols, function(x) ifelse(length(x) >
#											1, TRUE, FALSE))]
#		timevary <- as.character(c(timevary, unlist(detCols[toadd])))
#	}
#	dataLong <- reshape(data, timevary, direction = "long", sep = sep,
#			timevar = "timeINDEX", idvar = "idINDEX", ids = rownames(data))
#	## reorder dataLong to be in site-major order
#	dataLong <- dataLong[order(dataLong$idINDEX, dataLong$timeINDEX),]
##	mf <- model.frame(stateformula, dataLong)
#	browser()
##	y <- data[,yNames]
##	X <- model.matrix(stateformula, data)
#	obsCovs <- model.frame(detformula, dataLong)
#	siteCovs <- model.frame()
##	return(list(y=y, X=X, V=V, dataLong=dataLong))
#}




lambda2psi <- function(lambda)
{
if(any(lambda < 0))
	stop("lambda must be >= 0")
as.numeric(1 - exp(-lambda))
}


# Convert individual-level distance data to the transect-level format required by distsamp()

formatDistData <- function(distData, distCol, transectNameCol, dist.breaks)
{
    transects <- distData[,transectNameCol]
    M <- nlevels(transects)
    J <- length(dist.breaks) - 1
    y <- matrix(NA, M, J,
	   dimnames = list(levels(transects), paste("y", 1:J, sep=".")))
    for(i in 1:M) {
	   sub <- subset(distData, transects==rownames(y)[i])
	   y[i,] <- table(cut(sub[,distCol], dist.breaks, include.lowest=TRUE))
	   }
    return(y)
}



## Sight distance to perpendicular distance

sight2perpdist <- function(sightdist, sightangle)
{
    if(any(0 > sightangle | sightangle > 180))
	   stop("sightangle must be degrees in [0, 180]")
    sightdist * sin(sightangle * pi / 180)
}



SSE <- function(fit)
{
    sse <- sum(residuals(fit)^2, na.rm=TRUE)
    return(c(SSE=sse))
}



# For pcountOpen. Calculate time intervals acknowledging gaps due to NAs
# The first column indicates is time since first primary period + 1
formatDelta <- function(d, yna)
{
    M <- nrow(yna)
    T <- ncol(yna)
    d <- d - min(d, na.rm=TRUE) + 1
    dout <- matrix(NA, M, T)
    dout[,1] <- d[,1]
    dout[,2:T] <- t(apply(d, 1, diff))
    for(i in 1:M) {
        if(any(yna[i,]) & !all(yna[i,])) { # 2nd test for simulate
            last <- max(which(!yna[i,]))
            y.in <- yna[i, 1:last]
            d.in <- d[i, 1:last]
            if(any(y.in)) {
                for(j in last:2) { # first will always be time since 1
                    nextReal <- which(!yna[i, 1:(j-1)])
                    if(length(nextReal) > 0)
                        dout[i, j] <- d[i, j] - d[i, max(nextReal)]
                    else
                        dout[i, j] <- d[i, j] - 1
                    }
                }
            }
        }
    return(dout)
}












# Markov transition probs for pcountOpen
tranProbs <- function(Nr, omegaR, gammaR, deltaR, dynamicsR)
{
    if(any(Nr < 0))
        stop("N should be a non-negative integer")
    if(any(is.na(omegaR)) | any(is.na(gammaR)))
        stop("Missing values are not allowed in omega or gamma")
    if(length(omegaR) != 1 | length(gammaR) != 1 | length(deltaR) != 1)
        stop("omega, gamma, and delta must be scalars")
    if(any(is.na(deltaR)))
        stop("Delta cannot be NA")
    if(!dynamicsR %in% c("constant", "autoreg", "notrend"))
        stop("dynamics must be one of: constant, autoreg, or notrend")

    .Call("tranProbs",
        as.integer(Nr),
        as.double(omegaR),
        as.double(gammaR),
        as.integer(deltaR),
        as.character(dynamicsR),
        PACKAGE = "unmarked")
}


# The slow way
tranProbsR <- function(N, omega, gamma, delta, dynamics) {
    lN <- length(N)
    bpsum <- matrix(NA, lN, lN)
    for(j in 1:lN) {
        for(k in 1:lN) {
            cmin0 <- 0:min(N[j], N[k])
            gamma2 <- ifelse(identical(dynamics, "autoreg"), gamma*N[j], gamma)
            bpsum[k, j] <- sum(dbinom(cmin0, N[j], omega) *
                dpois(N[k]-cmin0, gamma2))
            }}
    bpsum <- t(t(bpsum) / colSums(bpsum))
    if(delta>1) {
        for(d in 2:delta) {
            bpsum <- bpsum %*% bpsum
            bpsum <- t(t(bpsum) / colSums(bpsum))
            }
        }
    return(bpsum)
    }




# much faster than apply(m, 1, prod) or exp(rowSums(log(m)))
rowProds <- function(m)
{
    .Call("rowProds",
        as.matrix(m),
        PACKAGE = "unmarked")
}




# ------------------------ fit statistics -------------------------------


