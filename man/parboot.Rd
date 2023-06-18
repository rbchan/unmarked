\name{parboot}
\alias{parboot}
\alias{parboot,unmarkedFitIDS-method}
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
  \item{seed}{set seed for reproducible bootstrap}
  \item{parallel}{logical (default = \code{TRUE}) indicating whether to compute 
    bootstrap on multiple cores, if present.  If \code{TRUE}, suppresses reporting
    of bootstrapped statistics.  Defaults to serial calculation when \code{nsim} < 100.
    Parallel computation is likely to be slower for simple models when \code{nsim} < ~500,
    but should speed up the bootstrap of more complicated models.}
  \item{ncores}{integer (default = one less than number of available cores) number of cores to
    use when bootstrapping in parallel.} 
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

\author{Richard Chandler \email{rbchan@uga.edu} and Adam Smith}

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
fitstats <- function(fm, na.rm=TRUE) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2, na.rm=na.rm)
    chisq <- sum((observed - expected)^2 / expected, na.rm=na.rm)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=na.rm)
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
