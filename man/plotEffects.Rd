\name{plotEffects}
\alias{plotEffects}
\alias{plotEffects-methods}
\alias{plotEffects,unmarkedFit-method}
\alias{plotEffectsData}
\alias{plotEffectsData-methods}
\alias{plotEffectsData,unmarkedFit-method}

\title{Plot marginal effects of covariates in unmarked models}

\description{This function generates a plot visualizing the effects
 of a single covariate on a parameter (e.g. occupancy, abundance) in an unmarked 
 model. If the covariate is numeric, the result is a line plot with an error
 ribbon where the x-axis is the range of the covariate and the y-axis is the 
 predicted parameter value. If the covariate is an R factor (i.e., categorical),
 the x-axis instead contains each unique value of the covariate. 

 All covariates in the model besides the one being plotted are held either at
 their median value (if they are numeric) or at their reference level (if they
 are factors).

 Some types of unmarked models may require additional arguments, which are passed
 to the matching \code{predict} method. For example, \code{unmarkedFitOccuMulti} 
 models require the \code{species} argument to be included in the function 
 call in order to work properly.

 If you want to customize a plot, the easiest approach is to get data 
 formatted for plotting using \code{plotEffectsData}, and use that. If
 you want to see and/or modify the code used by \code{plotEffects} to generate
 the default plots, run \code{getMethod("plotEffects", "unmarkedFit")} in
 the R console.
}

\usage{
\S4method{plotEffects}{unmarkedFit}(object, type, covariate, level=0.95, ...)
\S4method{plotEffectsData}{unmarkedFit}(object, type, covariate, level=0.95, ...)
}

\arguments{
  \item{object}{A fitted model inheriting class \code{unmarkedFit}}
  \item{type}{Submodel in which the covariate of interest can be found, for
    example \code{"state"} or \code{"det"}. This will depend on the fitted model}
  \item{covariate}{The name of the covariate to be plotted, as a character string}
  \item{level}{Confidence level for the error ribbons or bars}
  \item{...}{Other arguments passed to the \code{predict} function, required
    for some \code{unmarkedFit} types such as \code{unmarkedFitOccuMulti}}
} 

\value{A plot (\code{plotEffects} or a data frame (\code{plotEffectsData})
  containing values to be used in a plot.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\examples{

\dontrun{

# Simulate data and build an unmarked frame
set.seed(123)
dat_occ <- data.frame(x1=rnorm(500))
dat_p <- data.frame(x2=rnorm(500*5))

y <- matrix(NA, 500, 5)
z <- rep(NA, 500)

b <- c(0.4, -0.5, 0.3, 0.5)

re_fac <- factor(sample(letters[1:5], 500, replace=T))
dat_occ$group <- re_fac
re <- rnorm(5, 0, 1.2)
re_idx <- as.numeric(re_fac)

idx <- 1
for (i in 1:500){
  z[i] <- rbinom(1,1, plogis(b[1] + b[2]*dat_occ$x1[i] + re[re_idx[i]]))
  for (j in 1:5){
    y[i,j] <- z[i]*rbinom(1,1,
                    plogis(b[3] + b[4]*dat_p$x2[idx]))
    idx <- idx + 1
  }
}

umf <- unmarkedFrameOccu(y=y, siteCovs=dat_occ, obsCovs=dat_p)

# Fit model
(fm <- occu(~x2 ~x1 + group, umf))

# Plot marginal effects of various covariates
plotEffects(fm, "state", "x1")
plotEffects(fm, "state", "group")
plotEffects(fm, "det", "x2")

# Get raw data used for a plot
plotEffectsData(fm, "state", "group")

# See code used by plotEffects so you can edit it yourself and customize the plot
methods::getMethod("plotEffects", "unmarkedFit")
}
}
