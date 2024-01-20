\name{unmarkedEstimate-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedEstimate}
\alias{unmarkedEstimate-class}
\alias{summary,unmarkedEstimate-method}
\title{Class "unmarkedEstimate"}
\description{Contains parameter estimates, covariance matrix, and metadata}
\section{Objects from the Class}{
Creating these objects is done internally not by users. 
}
\section{Slots}{
	 \describe{
    \item{\code{name}:}{Object of class \code{"character"} storing parameter 
		names}
    \item{\code{short.name}:}{Object of class \code{"character"} storing 
		abbreviated parameter names}
    \item{\code{estimates}:}{Object of class \code{"numeric"}}
    \item{\code{covMat}:}{Object of class \code{"matrix"}}
    \item{\code{covMatBS}:}{Object of class \code{"matrix"}}
    \item{\code{fixed}:}{Object of class \code{"numeric"}}
    \item{\code{invlink}:}{Object of class \code{"character"}}
    \item{\code{invlinkGrad}:}{Object of class \code{"character"}}
    \item{\code{randomVarInfo}:}{Object of class \code{"list"}}
  }
}
\section{Methods}{
  \describe{
    \item{backTransform}{\code{signature(obj = "unmarkedEstimate")}}
    \item{coef}{\code{signature(object = "unmarkedEstimate")}}
    \item{confint}{\code{signature(object = "unmarkedEstimate")}}
    \item{linearComb}{\code{signature(obj = "unmarkedEstimate", 
		coefficients = "matrixOrVector")}}
    \item{SE}{\code{signature(obj = "unmarkedEstimate")}}
    \item{vcov}{\code{signature(object = "unmarkedEstimate")}}
	 }
}
\note{These methods are typically called within a call to a method for 
	\code{\link{unmarkedFit-class}}}
\examples{
showClass("unmarkedEstimate")
}
\keyword{classes}
