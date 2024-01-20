\name{unmarkedFrame-extractors}
\alias{getY}
\alias{siteCovs}
\alias{siteCovs<-}
\alias{obsCovs}
\alias{obsCovs<-}
\alias{yearlySiteCovs}
\alias{yearlySiteCovs<-}
\alias{obsToY}
\alias{obsToY<-}
\alias{getY,unmarkedFrame-method}
\alias{siteCovs,unmarkedFrame-method}
\alias{siteCovs<-,unmarkedFrame-method}
\alias{obsCovs,unmarkedFrame-method}
\alias{obsCovs<-,unmarkedFrame-method}
\alias{yearlySiteCovs,unmarkedMultFrame-method}
\alias{yearlySiteCovs<-,unmarkedMultFrame-method}
\alias{obsToY,unmarkedFrame-method}
\alias{obsToY<-,unmarkedFrame-method}

\title{Extract and replace methods for unmarkedFrame components}

\description{Methods for extracting and replacing components of \code{unmarkedFrame}s, 
  such as the observation data (\code{getY}), site covariates (\code{siteCovs}),
  and observation covariates (\code{obsCovs}).
}

\usage{
\S4method{getY}{unmarkedFrame}(object)
\S4method{siteCovs}{unmarkedFrame}(object)
\S4method{obsCovs}{unmarkedFrame}(object, matrices = FALSE)
\S4method{yearlySiteCovs}{unmarkedMultFrame}(object)
\S4method{obsToY}{unmarkedFrame}(object)
}

\arguments{
	\item{object}{An object of class \code{unmarkedFrame}}
  \item{matrices}{Logical; if FALSE covariates are returned as a long 
    data frame (MJ x ncovs). If true, a list of matrices (each M x J).}
}

\value{
  A matrix of observation data (\code{getY}), a data frame of covariates 
  (\code{siteCovs}, \code{obsCovs}, \code{yearlySiteCovs}), or the 
  \code{obsToY} matrix (\code{obsToY}).
}

\examples{

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
	obsCovs = mallard.obs)

# Extract the observations
getY(mallardUMF)

# Extract the observation covariates
obsCovs(mallardUMF)
obsCovs(mallardUMF, matrices = TRUE)

# Extract the site covariates
siteCovs(mallardUMF)

# Replace the existing site covariates
new_covs <- data.frame(x = rnorm(239))
siteCovs(mallardUMF) <- new_covs
siteCovs(mallardUMF)
}
