\name{predict-methods}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,ANY-method}
\alias{predict,unmarkedFit-method}
\alias{predict,unmarkedFitOccuFP-method}
\alias{predict,unmarkedFitOccuMulti-method}
\alias{predict,unmarkedFitOccuMS-method}
\alias{predict,unmarkedFitOccuTTD-method}
\alias{predict,unmarkedFitNmixTTD-method}
\alias{predict,unmarkedFitPCount-method}
\alias{predict,unmarkedFitColExt-method}
\alias{predict,unmarkedFitGMM-method}
\alias{predict,unmarkedFitGDS-method}
\alias{predict,unmarkedFitPCO-method}
\alias{predict,unmarkedFitDSO-method}
\alias{predict,unmarkedFitGDR-method}
\alias{predict,unmarkedFitIDS-method}
\alias{predict,unmarkedFitOccuRNMulti-method}
\alias{predict,unmarkedFitList-method}
\alias{predict,unmarkedRanef-method}
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
\item{\code{signature(object = "unmarkedRanef")}}{
Use this method to generate the empirical Bayes posterior predictive distribution
for functions of the random variables (latent abundance or occurrence).

In addition to the output object from \code{ranef}, you must also supply a 
custom function to argument \code{func}. The function must take as input a matrix
with dimensions M x T, where M is the number of sites and T is the number of 
primary periods (T=1 for single-season models). The output of this function should 
be a vector or matrix containing the derived parameters of interest. 

You may also manually set the number of draws from the posterior predictive
distribution with argument \code{nsims}; the default is 100.

The output of \code{predict} will be a vector or array with one more dimension
than the output of the function supplied \code{func}, corresponding to the number
of draws requested \code{nsims}. For example, if \code{func}
outputs a scalar, the output of \code{predict} will be a vector with length
equal to \code{nsims}. If \code{func} outputs a 3x2 matrix, the output of
\code{predict} will be an array with dimensions 3x2x\code{nsims}.
See \code{\link{ranef}} for an example.

Alternatively, you can use the \code{\link{posteriorSamples}} function on the
\code{ranef} output object to obtain the full posterior predictive distribution. 
This is useful if you are having trouble designing your custom function or if 
you want to obtain multiple different derived parameters from the same posterior
predictive distribution.

}
}}
\keyword{methods}
