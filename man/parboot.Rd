\name{parboot}
\alias{parboot}
\alias{plot,parboot,missing-method}
\alias{show,parboot-method}
\title{Parametric bootstrap method for fitted models inheriting class.}
\description{Simulate datasets from a fitted model, refit the model, and
  generate a sampling distribution for a user-specified fit-statistic.}
\arguments{
  \item{object}{a fitted model inheriting class "unmarkedFit"}
  \item{statistic}{a function returning a vector of fit-statistics.
    First argument must be the fitted model.
    Default is sum of squared residuals.}
  \item{nsim}{number of bootstrap replicates}
  \item{report}{print fit statistic every 'report' iterations during resampling}
  \item{...}{Additional arguments to be passed to statistic}}
\details{This function simulates datasets based upon a fitted model,
  refits the model, and evaluates a user-specified fit-statistic for each
  simulation. Comparing this sampling distribution to the observed statistic
  provides a means of evaluating goodness-of-fit or assessing uncertainty in
  a quantity of interest.}
\value{
  An object of class parboot with three slots:
  \item{call}{parboot call}
  \item{t0}{Numeric vector of statistics for original fitted model.}
  \item{t.star}{nsim by length(t0) matrix of statistics for each simulation fit.}}

\author{Richard Chandler \email{rbchan@uga.edu}}

\seealso{
  \code{\link{ranef}}
}

\examples{

data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20))
lengths <- linetran$Length

ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths*1000, survey = "line", unitsIn = "m")
    })

# Fit a model
(fm <- distsamp(~area ~habitat, ltUMF))

# Function returning three fit-statistics.
fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2)
    chisq <- sum((observed - expected)^2 / expected)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
    }

(pb <- parboot(fm, fitstats, nsim=25, report=1))
plot(pb, main="")


# Finite-sample inference for a derived parameter.
# Population size in sampled area

Nhat <- function(fm) {
    sum(bup(ranef(fm, K=50)))
    }

set.seed(345)
(pb.N <- parboot(fm, Nhat, nsim=25, report=5))

# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(fm, K=50)))



}