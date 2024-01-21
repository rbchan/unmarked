\name{unmarkedFrameOccuFP}

\title{Organize data for the single season occupancy models fit by occuFP}

\alias{unmarkedFrameOccuFP}

\usage{unmarkedFrameOccuFP(y, siteCovs=NULL, obsCovs=NULL, type, mapInfo)}

\description{Organizes detection, non-detection data along with the covariates.
    This S4 class is required by the data argument of \code{\link{occu}} and
    \code{\link{occuRN}}}

\arguments{
    \item{y}{An RxJ matrix of the detection, non-detection data, where R is the
        number of sites, J is the maximum number of sampling periods per site.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of
        covariates that vary within sites, or a \code{\link{data.frame}} with
        RxJ rows in site-major order.}
    \item{type}{A vector with 3 values designating the number of occassions
	where data is of type 1, type 2, and type 3 - see \code{\link{occuFP}}
	for more details about data types.}
    \item{mapInfo}{Currently ignored}

}

\details{
    unmarkedFrameOccuFP is the S4 class that holds data to be passed
    to the \code{\link{occu}} and \code{\link{occuRN}} model-fitting function.
}

\value{an object of class unmarkedFrameOccuFP}

\seealso{\code{\link{unmarkedFrame-methods}}, \code{\link{occuFP}}}

\examples{

n = 100
o = 10
o1 = 5
y = matrix(0,n,o)
p = .7
r = .5
fp = 0.05
y[1:(n*.5),(o-o1+1):o] <- rbinom((n*o1*.5),1,p)
y[1:(n*.5),1:(o-o1)] <- rbinom((o-o1)*n*.5,1,r)
y[(n*.5+1):n,(o-o1+1):o] <- rbinom((n*o1*.5),1,fp)
type <- c((o-o1),o1,0)  ### vector with the number of each data type
site <- c(rep(1,n*.5*.8),rep(0,n*.5*.2),rep(1,n*.5*.2),rep(0,n*.8*.5))
occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
site <- data.frame(habitat = site)
occ <- list(METH = occ)

umf1 <- unmarkedFrameOccuFP(y,site,occ, type = type)

m1 <- occuFP(detformula = ~ METH, FPformula = ~1, stateformula = ~ habitat, data = umf1)

}
