\name{predict-methods}
\docType{methods}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,unmarkedFit-method}
\alias{predict,unmarkedFitPCount-method}
\alias{predict,unmarkedFitColExt-method}
\alias{predict,unmarkedFitGMM-method}
\alias{predict,unmarkedFitGDS-method}
\alias{predict,unmarkedFitPCO-method}
\alias{predict,unmarkedFitList-method}
\title{ Methods for Function predict in Package `unmarked' }
\description{
These methods return predicted values from fitted model objects.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFit")}}{
"type" must be either `state' or `det'.
}
\item{\code{signature(object = "unmarkedFitColExt")}}{
"type" must be 'psi', 'col', 'ext', or 'det'.
}
\item{\code{signature(object = "unmarkedFitGMM")}}{
"type" must be 'lambda', 'psi', 'det'
}
\item{\code{signature(object = "unmarkedFitList")}}{
"type" depends upon the fitted models
}
}}
\keyword{methods}
