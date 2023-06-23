\name{crossVal}
\alias{crossVal}
\alias{crossVal-methods}
\alias{crossVal,unmarkedFit-method}
\alias{crossVal,unmarkedFitList-method}
\alias{unmarkedCrossVal-class}
\alias{unmarkedCrossValList-class}
\alias{show,unmarkedCrossVal-method}
\alias{show,unmarkedCrossValList-method}

\title{Cross-validation methods for fitted unmarked models and fit lists}

\description{Test predictive accuracy of fitted models using several 
  cross-validation approaches. The dataset is divided by site only into folds or
  testing and training datasets (i.e., encounter histories within sites are
  never split up).}

\usage{
\S4method{crossVal}{unmarkedFit}(
  object, method=c("Kfold","holdout","leaveOneOut"),
  folds=10, holdoutPct=0.25, statistic=RMSE_MAE, parallel=FALSE, ncores, ...)
\S4method{crossVal}{unmarkedFitList}(
  object, method=c("Kfold","holdout","leaveOneOut"),
  folds=10, holdoutPct=0.25, statistic=RMSE_MAE, parallel=FALSE, ncores, 
  sort = c("none", "increasing", "decreasing"), ...)
}

\arguments{
  \item{object}{A fitted model inheriting class \code{unmarkedFit} or a list of 
    fitted models with class \code{unmarkedFitList}}
  \item{method}{Cross validation method to use as string. Valid options are 
    \code{"Kfold"}, \code{"holdout"}, or \code{"leaveOneOut"}}
  \item{folds}{Number of folds to use for k-fold cross validation}
  \item{holdoutPct}{Proportion of dataset (value between 0-1) to use as the 
    "holdout" or "test" set, for the holdout method}
  \item{statistic}{Function that calculates statistics for each fold. The
    function must take an \code{unmarkedFit} object as the first argument and
    return a named numeric vector with statistic value(s). The default function
    \code{RMSE_MAE} returns root-mean-square error and mean absolute error.
    See \code{unmarked:::RMSE_MAE} for an example of correct statistic
    function structure.}
  \item{parallel}{If \code{TRUE}, run folds in parallel. This may speed up 
    cross-validation if the unmarked model takes a long time to fit or you have 
    a large number of sites and are using leave-one-out cross-validation.}
  \item{ncores}{Number of parallel cores to use.}
  \item{sort}{If doing cross-validation on a \code{fitList}, you can optionally
    sort the resulting table(s) of statistic values for each model.}
  \item{...}{Other arguments passed to the statistic function.}
} 

\value{\code{unmarkedCrossVal} or \code{unmarkedCrossValList} object containing
  calculated statistic values for each fold.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{fitList}},
  \code{\link{unmarkedFit}}
}

\examples{

\dontrun{
#Get data
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))    
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

#Fit occupancy model
fm <- occu(~ obsvar1 ~ 1, pferUMF)

#k-fold cross validation with 10 folds
(kfold = crossVal(fm, method="Kfold", folds=10))

#holdout method with 25% of data in test set and 75% in train set
(holdout = crossVal(fm,method='holdout', holdoutPct=0.25))

#Leave-one-out method
(leave = crossVal(fm, method='leaveOneOut'))

#Fit a second model and combine into a fitList
fm2 <- occu(~1 ~1, pferUMF)
fl <- fitList(fm2,fm)

#Cross-validation for all fits in fitList using holdout method
(cvlist <- crossVal(fl, method='holdout'))

}
}
