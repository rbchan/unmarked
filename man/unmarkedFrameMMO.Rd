\name{unmarkedFrameMMO}
\title{Create an object of class unmarkedFrameMMO that contains data used by multmixOpen.}
\alias{unmarkedFrameMMO}
\usage{unmarkedFrameMMO(y, siteCovs=NULL, obsCovs=NULL, yearlySiteCovs=NULL, 
       numPrimary, type, primaryPeriod)}

\description{Organizes count data and experimental design information 
  from multiple primary periods along with associated covariates. This S4 class 
  is required by the data argument of \code{\link{multmixOpen}}}

\arguments{
    \item{y}{An MxJT matrix of the repeated count data, where M is the number
        of sites (i.e., points or transects), J is the number of distance classes
        and T is the maximum number of primary sampling periods per site}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of covariates
        that vary within sites, or a \code{\link{data.frame}} with MxJT rows in
        site-major order.}
    \item{yearlySiteCovs}{Either a named list of MxT \code{\link{data.frame}}s,
        or a site-major \code{\link{data.frame}} with MT rows and 1 column per
        covariate}
    \item{numPrimary}{Maximum number of observed primary periods for each site}
    \item{type}{Either "removal" for removal sampling, "double" for standard 
        double observer sampling, or "depDouble" for dependent double observer 
        sampling}
    \item{primaryPeriod}{An MxJT matrix of integers indicating the primary period 
      of each observation}
}
\details{
  \code{unmarkedFrameMMO} is the S4 class that holds data to be passed
  to the \code{\link{multmixOpen}} model-fitting function.

  Options for the detection process (\code{type}) include equal-interval removal 
  sampling (\code{"removal"}), double observer sampling (\code{"double"}), or 
  dependent double-observer sampling (\code{"depDouble"}). Note
  that unlike the related functions \code{\link{multinomPois}} and
  \code{\link{gmultmix}}, custom functions for the detection process (i.e., 
  \code{piFun}s) are not supported. To request additional options contact the author.

  When gamma or omega are modeled using year-specific covariates, the
  covariate data for the final year will be ignored; however,
  they must be supplied.

  If the time gap between primary periods is not constant, an M by T
  matrix of integers should be supplied using the \code{primaryPeriod} argument.  
}

\value{an object of class unmarkedFrameMMO}
\seealso{\code{\link{unmarkedFrame-methods}}, \code{\link{multmixOpen}}}

\examples{

  #Generate some data 
  set.seed(123)
  lambda=4; gamma=0.5; omega=0.8; p=0.5
  M <- 100; T <- 5
  y <- array(NA, c(M, 3, T))
  N <- matrix(NA, M, T)
  S <- G <- matrix(NA, M, T-1)

  for(i in 1:M) {
    N[i,1] <- rpois(1, lambda)
    y[i,1,1] <- rbinom(1, N[i,1], p)    # Observe some
    Nleft1 <- N[i,1] - y[i,1,1]         # Remove them
    y[i,2,1] <- rbinom(1, Nleft1, p)   # ...
    Nleft2 <- Nleft1 - y[i,2,1]
    y[i,3,1] <- rbinom(1, Nleft2, p)

    for(t in 1:(T-1)) {
      S[i,t] <- rbinom(1, N[i,t], omega)
      G[i,t] <- rpois(1, gamma)
      N[i,t+1] <- S[i,t] + G[i,t]
      y[i,1,t+1] <- rbinom(1, N[i,t+1], p)    # Observe some
      Nleft1 <- N[i,t+1] - y[i,1,t+1]         # Remove them
      y[i,2,t+1] <- rbinom(1, Nleft1, p)   # ...
      Nleft2 <- Nleft1 - y[i,2,t+1]
      y[i,3,t+1] <- rbinom(1, Nleft2, p)
    }
  }
  y=matrix(y, M)
  
  #Create some random covariate data
  sc <- data.frame(x1=rnorm(100))

  #Create unmarked frame
  umf <- unmarkedFrameMMO(y=y, numPrimary=5, siteCovs=sc, type="removal")
  
  summary(umf)

}
