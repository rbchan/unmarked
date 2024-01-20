\name{unmarkedPower-methods}
\alias{unmarkedPower-methods}
\alias{unmarkedPower-class}
\alias{summary,unmarkedPower-method}
\alias{update,unmarkedPower-method}

\title{Methods for unmarkedPower objects}

\description{Various functions to summarize and update unmarkedPower objects}

\usage{
\S4method{summary}{unmarkedPower}(object, ...)
\S4method{update}{unmarkedPower}(object, ...)
}

\arguments{
  \item{object}{An object of class \code{unmarkedPower} created with the 
    \code{powerAnalysis} function}
  \item{...}{For \code{update}, arguments to change in the updated power analysis.
    Not used by \code{summary}}
}

\value{
  For \code{summary}, summary output is printed to the console.
  For \code{update}, a new \code{powerAnalysis} object corresponding to the 
  new arguments provided.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{powerAnalysis}}
}

\examples{

\dontrun{

# Simulate an occupancy dataset
forms <- list(state=~elev, det=~1)
coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
design <- list(M=300, J=8) # 300 sites, 8 occasions per site
occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)

# Fit occupancy model to simulated data
template_model <- occu(~1~elev, occu_umf)

# Set desired effect sizes to pass to coefs
effect_sizes <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

# Run power analysis
pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05)

# Look at summary
summary(pa)

# Update the analysis with new arguments
(pa2 <- update(pa, alpha=0.01))

}
}
