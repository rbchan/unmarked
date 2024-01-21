\name{unmarkedFrameOccuMulti}

\title{Organize data for the multispecies occupancy model fit by occuMulti}

\alias{unmarkedFrameOccuMulti}

\usage{unmarkedFrameOccuMulti(y, siteCovs=NULL, obsCovs=NULL, 
                              maxOrder, mapInfo)}

\description{Organizes detection, non-detection data for multiple species along 
    with the covariates. This S4 class is required by the data argument 
    of \code{\link{occuMulti}}}

\arguments{
    \item{y}{A list (optionally a named list) of length S where each element 
        is an MxJ matrix of the detection, non-detection data for one species, 
        where M is the number of sites, J is the maximum number of sampling 
        periods per site, and S is the number of species in the analysis.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of 
        covariates that vary within sites, or a \code{\link{data.frame}} with 
        MxJ rows in site-major order.}
    \item{maxOrder}{Optional; specify maximum interaction order. Defaults to 
        number of species (all possible interactions). Reducing this value may
        speed up creation of unmarked frame if you aren't interested in 
        higher-order interactions.}
    \item{mapInfo}{Currently ignored}
  
}

\details{
    unmarkedFrameOccuMulti is the S4 class that holds data to be passed 
    to the \code{\link{occuMulti}} model-fitting function.
}

\value{an object of class unmarkedFrameOccuMulti}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarkedFrame-methods}}, \code{\link{occuMulti}}}

\examples{

# Fake data
S <- 3 # number of species
M <- 4 # number of sites
J <- 3 # number of visits

y <- list(matrix(rbinom(M*J,1,0.5),M,J), # species 1
          matrix(rbinom(M*J,1,0.5),M,J), # species 2
          matrix(rbinom(M*J,1,0.2),M,J)) # species 3

site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

umf <- unmarkedFrameOccuMulti(y=y, siteCovs=site.covs, 
    obsCovs=NULL)   # organize data
umf                     # look at data
summary(umf)            # summarize      
plot(umf)               # visualize
#fm <- occu(~1 ~1, umf)  # fit a model


}
