\name{unmarkedFramePCO}
\title{Create an object of class unmarkedFramePCO that contains data used by pcountOpen.}
\alias{unmarkedFramePCO}
\usage{unmarkedFramePCO(y, siteCovs=NULL, obsCovs=NULL, yearlySiteCovs, mapInfo,
    numPrimary, primaryPeriod)}
\description{Organizes repeated count data along with the covariates and possibly the dates on which each survey was conducted. This S4 class is required by the data argument of \code{\link{pcountOpen}}}
\arguments{
    \item{y}{An MxJT matrix of the repeated count data, where M is the number
        of sites, J is the maximum number of secondary sampling periods per site
        and T is the maximum number of primary sampling periods per site.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of covariates
        that vary within sites, or a \code{\link{data.frame}} with MxJT rows in
        site-major order.}
    \item{yearlySiteCovs}{Either a named list of MxT \code{\link{data.frame}}s,
        or a site-major \code{\link{data.frame}} with MT rows and 1 column per
        covariate.}
    \item{mapInfo}{Currently ignored}
    \item{numPrimary}{Maximum number of observed primary periods for each site}
    \item{primaryPeriod}{matrix of integers indicating the primary period of each survey.}

}
\details{unmarkedFramePCO is the S4 class that holds data to be passed
to the \code{\link{pcountOpen}} model-fitting function.

The unmarkedFramePCO class is similar to the unmarkedFramePCount class except
that it contains the dates for each survey, which needs to be supplied .

}

\value{an object of class unmarkedFramePCO}
\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}},
    \code{\link{pcountOpen}}}
\examples{

# Repeated count data with 5 primary periods and
# no secondary sampling periods (ie J==1)
y1 <- matrix(c(
    0, 2, 3, 2, 0,
    2, 2, 3, 1, 1,
    1, 1, 0, 0, 3,
    0, 0, 0, 0, 0), nrow=4, ncol=5, byrow=TRUE)

# Site-specific covariates
sc1 <- data.frame(x1 = 1:4, x2 = c('A','A','B','B'))

# Observation-specific covariates
oc1 <- list(
    x3 = matrix(1:5, nrow=4, ncol=5, byrow=TRUE),
    x4 = matrix(letters[1:5], nrow=4, ncol=5, byrow=TRUE))

# Primary periods of surveys
primaryPeriod1 <- matrix(as.integer(c(
    1, 2, 5, 7, 8,
    1, 2, 3, 4, 5,
    1, 2, 4, 5, 6,
    1, 3, 5, 6, 7)), nrow=4, ncol=5, byrow=TRUE)


# Create the unmarkedFrame
umf1 <- unmarkedFramePCO(y=y1, siteCovs=sc1, obsCovs=oc1, numPrimary=5,
    primaryPeriod=primaryPeriod1)

# Take a look
umf1
summary(umf1)






# Repeated count data with 4 primary periods and
# no 2 secondary sampling periods (ie J=2)
y2 <- matrix(c(
    0,0,  2,2,  3,2,  2,2,
    2,2,  2,1,  3,2,  1,1,
    1,0,  1,1,  0,0,  0,0,
    0,0,  0,0,  0,0,  0,0), nrow=4, ncol=8, byrow=TRUE)


# Site-specific covariates
sc2 <- data.frame(x1 = 1:4, x2 = c('A','A','B','B'))

# Observation-specific covariates
oc2 <- list(
    x3 = matrix(1:8, nrow=4, ncol=8, byrow=TRUE),
    x4 = matrix(letters[1:8], nrow=4, ncol=8, byrow=TRUE))

# Yearly-site covariates
ysc <- list(
    x5 = matrix(c(
        1,2,3,4,
        1,2,3,4,
        1,2,3,4,
        1,2,3,4), nrow=4, ncol=4, byrow=TRUE))

# Primary periods of surveys
primaryPeriod2 <- matrix(as.integer(c(
    1,2,5,7,
    1,2,3,4,
    1,2,4,5,
    1,3,5,6)), nrow=4, ncol=4, byrow=TRUE)

# Create the unmarkedFrame
umf2 <- unmarkedFramePCO(y=y2, siteCovs=sc2, obsCovs=oc2,
    yearlySiteCovs=ysc,
    numPrimary=4, primaryPeriod=primaryPeriod2)

# Take a look
umf2
summary(umf2)


}
