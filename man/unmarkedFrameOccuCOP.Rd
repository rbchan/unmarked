\name{unmarkedFrameOccuCOP}
\alias{unmarkedFrameOccuCOP}
\alias{getL}
\alias{getL,unmarkedFrameOccuCOP-method}

\title{Organize data for the occupancy model using count data fit by \code{occuCOP}}

\usage{unmarkedFrameOccuCOP(y, L, siteCovs = NULL, obsCovs = NULL)}

\description{Organizes count data along with the covariates. The \linkS4class{unmarkedFrame} S4 class required by the \code{data} argument of \code{\link{occuCOP}}.}

\arguments{
  
  \item{y}{An MxJ matrix of the count data, where M is the number of sites, J is the maximum number of observation periods (sampling occasions, transects, discretised sessions...) per site.}

  \item{L}{An MxJ matrix of the length of the observation periods. For example, duration of the sampling occasion in hours, duration of the discretised session in days, or length of the transect in meters.}
  
  \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the site level. This should have M rows and one column per covariate}
  
  \item{obsCovs}{A named list of dataframes of dimension MxJ, with one dataframe per covariate that varies between sites and observation periods}
  
  %% \item{mapInfo}{Currently ignored}
}

\details{
    unmarkedFrameOccuCOP is the \linkS4class{unmarkedFrame} S4 class that holds data to be passed to the \code{\link{occuCOP}} model-fitting function.
}

\value{an object of class \code{unmarkedFrameOccuCOP}}

\seealso{
  \code{\link{unmarkedFrame-class}}, 
  \code{\link{unmarkedFrame}}, 
  \code{\link{occuCOP}}
}

\examples{
# Fake data
M <- 4 # Number of sites
J <- 3 # Number of observation periods

# Count data
(y <- matrix(
  c(1, 3, 0,
    0, 0, 0,
    2, 0, 5,
    1, NA, 0),
  nrow = M,
  ncol = J,
  byrow = TRUE
))

# Length of observation periods
(L <- matrix(
  c(1, 3, NA,
    2, 2, 2,
    1, 2, 1,
    7, 1, 3),
  nrow = M,
  ncol = J,
  byrow = TRUE
))

# Site covariates
(site.covs <- data.frame(
  "elev" = rexp(4),
  "habitat" = factor(c("forest", "forest", "grassland", "grassland"))
))

# Observation covariates (as a list)
(obs.covs.list <- list(
  "rain" = matrix(rexp(M * J), nrow = M, ncol = J),
  "wind" = matrix(
    sample(letters[1:3], replace = TRUE, size = M * J),
    nrow = M, ncol = J)
))

# Organise data in a unmarkedFrameOccuCOP object
umf <- unmarkedFrameOccuCOP(
  y = y,
  L = L,
  siteCovs = site.covs,
  obsCovs = obs.covs.list
)

# Extract L
getL(umf)

# Look at data
print(umf) # Print the whole data set
print(umf[1, 2]) # Print the data of the 1st site, 2nd observation
summary(umf) # Summarise the data set
plot(umf) # Plot the count of detection events


# L is optional, if absent, it will be replaced by a MxJ matrix of 1
unmarkedFrameOccuCOP(
  y = y,
  siteCovs = site.covs,
  obsCovs = obs.covs.list
)

# Covariates are optional
unmarkedFrameOccuCOP(y = y)
}
