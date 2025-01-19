\name{linearComb-methods}
\docType{methods}
\alias{linearComb}
\alias{linearComb-methods}
\alias{linearComb,unmarkedEstimate,matrixOrVector-method}
\alias{linearComb,unmarkedFit,matrixOrVector-method}
\alias{linearComb,unmarkedFitOccuRNMulti,matrixOrVector-method}
\alias{show,unmarkedLinComb-method}
\title{Methods for Function linearComb in Package `unmarked'}
\description{
Methods for function \code{linearComb} in Package `unmarked'
}
\section{Methods}{
\describe{

\item{obj = "unmarkedEstimate", coefficients = "matrixOrVector"}{Typically 
	called internally}

\item{obj = "unmarkedFit", coefficients = "matrixOrVector"}{Returns 
  linear combinations of parameters from a fitted model.  Coefficients
  are supplied through coefficients.  The required argument \code{type}
  specifies which model estimate to use.  You can use
  \code{names(fittedmodel)} to view possible values for the \code{type} argument.}
}}
\keyword{methods}
\examples{
data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])), type = "removal")
fm <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)
linearComb(fm, c(1, 0.5, 0.5), type = "state")
linearComb(fm, matrix(c(1, 0.5, 0.5, 1, 0, 0, 1, 0, 0.5), 3, 3,
  byrow=TRUE), type="state")
}
