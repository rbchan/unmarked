\name{unmarkedRanef-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedRanef-class}
\alias{bup}
\alias{bup,unmarkedRanef-method}
\alias{confint,unmarkedRanef-method}
\alias{plot,unmarkedRanef,missing-method}
\alias{show,unmarkedRanef-method}

\title{Class \code{"unmarkedRanef"}}

\description{
Stores the estimated posterior distributions of
the latent abundance or occurrence variables.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{ranef}}.
}

\section{Slots}{
  \describe{
    \item{\code{post}:}{An \code{\link{array}} with nSites rows and Nmax
  (K+1) columns and nPrimaryPeriod slices }
  }
}

\section{Methods}{
  \describe{
    \item{bup}{\code{signature(object = "unmarkedRanef")}: Extract the
    Best Unbiased Predictors (BUPs) of the latent variables (abundance
  or occurrence state). Either the posterior mean or median can be
  requested using the \code{stat} argument.}
    \item{confint}{\code{signature(object = "unmarkedRanef")}: Compute
      confidence intervals.}
    \item{plot}{\code{signature(x = "unmarkedRanef", y = "missing")}:
      Plot the posteriors using \code{\link{xyplot}}}
    \item{show}{\code{signature(object = "unmarkedRanef")}: Display the
    modes and confidence intervals}
}
}

\section{Warnings}{
  Empirical Bayes methods can underestimate the variance of the
  posterior distribution because they do not account for uncertainty in
  the hyperparameters (lambda or psi). However, we have not found this
  to be the case in a set of limited simulations studies that can be
  found in the inst/unitTests directory of unmarked. Simulations do,
  however, indicate that the posterior mode can exhibit slight (3-5
  percent) negatively bias as a point
  estimator of site-specific abundance. It appears to be safer to use
  the posterior mean even though this will not be an integer in general.
}

\references{

  Laird, N.M. and T.A. Louis. 1987. Empirical Bayes confidence intervals
  based on bootstrap samples. Journal of the American Statistical
  Association 82:739--750.

  Carlin, B.P and T.A Louis. 1996. Bayes and Empirical Bayes Methods for
  Data Analysis. Chapman and Hall/CRC.

  Royle, J.A and R.M. Dorazio. 2008. Hierarchical Modeling and Inference
  in Ecology. Academic Press.

  }


\seealso{
\code{\link{ranef}}
}
\examples{
showClass("unmarkedRanef")
}
\keyword{classes}
