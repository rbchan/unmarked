\name{unmarkedPowerList}
\alias{unmarkedPowerList}
\alias{unmarkedPowerList,list-method}
\alias{unmarkedPowerList,unmarkedFit-method}
\alias{unmarkedPowerList-class}
\alias{unmarkedPowerList-methods}
\alias{show,unmarkedPowerList-method}
\alias{summary,unmarkedPowerList-method}
\alias{plot,unmarkedPowerList,ANY-method}

\title{Create or summarize a series of unmarked power analyses}

\description{
  A list of power analyses created with \code{powerAnalysis} can be combined 
  using \code{unmarkedPowerList}, allowing comparison e.g. between different
  study designs/sample sizes. Additionally an \code{unmarkedPowerList} can be
  created directly from an \code{unmarkedFit} template model by specifying
  a series of study designs (number of sites, number of observations) 
  as a \code{data.frame}. A series of methods for \code{unmarkedPowerList}
  objects are available including a \code{plot} method.
}

\usage{
\S4method{unmarkedPowerList}{list}(object, ...)
\S4method{unmarkedPowerList}{unmarkedFit}(object, coefs, design, alpha=0.05,
                                          nulls=list(), nsim=100, parallel=FALSE, ...)
\S4method{show}{unmarkedPowerList}(object)
\S4method{summary}{unmarkedPowerList}(object, ...)
\S4method{plot}{unmarkedPowerList,ANY}(x, power=NULL, param=NULL, ...)
}

\arguments{
  \item{object,x}{A \code{list} of \code{unmarkedPower} objects, a fitted model 
    inheriting class \code{unmarkedFit}, or an \code{unmarkedPowerList} object,
    depending on the method
  }
  \item{coefs}{A named list of effect sizes, see documentation for
    \code{powerAnalysis}}
  \item{design}{A \code{data.frame} with one row per study design to test, and
    at least 2 named columns: \code{M} for number of sites and \code{J} for
    number of observations. If you have >1 primary period a \code{T} column
    must also be provided}
  \item{alpha}{Type I error rate}
  \item{nulls}{If provided, a list matching the structure of \code{coefs} which
    defines the null hypothesis value for each parameter. By default the null
    is 0 for all parameters.
  }
  \item{nsim}{The number of simulations to run for each scenario/study design}
  \item{parallel}{If \code{TRUE}, run simulations in parallel}
  \item{power}{When plotting, the target power. Draws a horizontal line
    at a given value of power on the plot}
  \item{param}{When plotting, the model parameter to plot power vs. sample size for.
    By default this is the first parameter (which is usually an intercept,
    so not very interesting)}
  \item{...}{Not used}
} 

\value{A \code{unmarkedPowerList} object, a summary of the object in the console,
      or a summary plot, depending on the method}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{powerAnalysis}}
}

\examples{

\dontrun{

# Simulate an occupancy dataset and build template model
forms <- list(state=~elev, det=~1)
coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
design <- list(M=300, J=8) # 300 sites, 8 occasions per site
occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)
template_model <- occu(~1~elev, occu_umf)

# Generate two power analysis
effect_sizes <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05)
pa2 <- powerAnalysis(template_model, effect_sizes, design=list(M=100,J=2))

# Build unmarkedPowerList and look at summary
(pl <- unmarkedPowerList(list(pa,pa2)))

# Run a bunch of power analyses for different scenarios all at once
scenarios <- expand.grid(M=c(50,200,400),
                         J=c(3,5,8))
(pl2 <- unmarkedPowerList(template_model, effect_sizes, design=scenarios, nsim=20))

# Look at summary plot for elev effect
plot(pl2, power=0.8, param='elev')

}
}
