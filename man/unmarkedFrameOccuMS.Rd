\name{unmarkedFrameOccuMS}

\title{Organize data for the multi-state occupancy model fit by occuMS}

\alias{unmarkedFrameOccuMS}

\usage{unmarkedFrameOccuMS(y, siteCovs=NULL, obsCovs=NULL, 
                           numPrimary=1, yearlySiteCovs=NULL)}

\description{Organizes multi-state occupancy data (currently single-season only) 
  along with covariates. This S4 class is required by the data argument 
    of \code{\link{occuMS}}}

\arguments{
    \item{y}{An MxR matrix of multi-state occupancy data for a species, 
        where M is the number of sites and R is the maximum number of 
        observations per site (across all primary and secondary periods, if 
        you have multi-season data). Values in \code{y} should be integers ranging from
        0 (non-detection) to the number of total states - 1. For example, if you
        have 3 occupancy states, \code{y} should contain only values 0, 1, or 2.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of 
        covariates that vary within sites, or a \code{\link{data.frame}} with 
        MxR rows in the ordered by site-observation (if single-season) or 
        site-primary period-observation (if multi-season).}
    \item{numPrimary}{Number of primary time periods (e.g. seasons) for the 
        dynamic or multi-season version of the model. There should be
        an equal number of secondary periods in each primary period.}
    \item{yearlySiteCovs}{A data frame with one column per covariate that varies 
        among sites and primary periods (e.g. years). It should have MxT rows
        where M is the number of sites and T the number of primary periods,
        ordered by site-primary period. These covariates only used for dynamic 
        (multi-season) models.}
  
}

\details{
    unmarkedFrameOccuMS is the S4 class that holds data to be passed 
    to the \code{\link{occuMS}} model-fitting function.
}

\value{an object of class unmarkedFrameOccuMS}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarkedFrame-methods}}, \code{\link{occuMS}}}

\examples{

# Fake data
#Parameters
N <- 100; J <- 3; S <- 3
psi <- c(0.5,0.3,0.2)
p11 <- 0.4; p12 <- 0.25; p22 <- 0.3

#Simulate state
z <- sample(0:2, N, replace=TRUE, prob=psi)

#Simulate detection
y <- matrix(0,nrow=N,ncol=J)
for (n in 1:N){
  probs <- switch(z[n]+1,
                  c(0,0,0),
                  c(1-p11,p11,0),
                  c(1-p12-p22,p12,p22))
  
  if(z[n]>0){
    y[n,] <- sample(0:2, J, replace=TRUE, probs)
  }
}

#Covariates
site_covs <- as.data.frame(matrix(rnorm(N*2),ncol=2)) # nrow = # of sites
obs_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2)) # nrow = N*J

#Build unmarked frame
umf <- unmarkedFrameOccuMS(y=y,siteCovs=site_covs,obsCovs=obs_covs)

umf                     # look at data
summary(umf)            # summarize      
plot(umf)               # visualize
umf@numStates           # check number of occupancy states detected
}
