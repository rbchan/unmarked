\name{getP-methods}
\docType{methods}
\alias{getP}
\alias{getP-methods}
\alias{getP,unmarkedFit-method}
\alias{getP,unmarkedFitOccuFP-method}
\alias{getP,unmarkedFitOccuMulti-method}
\alias{getP,unmarkedFitOccuMS-method}
\alias{getP,unmarkedFitOccuTTD-method}
\alias{getP,unmarkedFitDS-method}
\alias{getP,unmarkedFitMPois-method}
\alias{getP,unmarkedFitColExt-method}
\alias{getP,unmarkedFitPCO-method}
\alias{getP,unmarkedFitGMM-method}
\alias{getP,unmarkedFitGDS-method}
\alias{getP,unmarkedFitGPC-method}
\alias{getP,unmarkedFitDSO-method}
\alias{getP,unmarkedFitMMO-method}
\alias{getP,unmarkedFitGDR-method}
\alias{getP,unmarkedFitIDS-method}
\alias{getP,unmarkedFitGOccu-method}
\alias{getP,unmarkedFitOccuCOP-method}
\alias{getP,unmarkedFitOccuRNMulti-method}

\title{Methods for Function getP in Package `unmarked'}
\description{
Methods for function \code{getP} in Package `unmarked'. These methods return a matrix of the back-transformed detection parameter (\eqn{p} the detection probability or \eqn{\lambda} the detection rate, depending on the model). The matrix is of dimension MxJ, with M the number of sites and J the number of sampling periods; or of dimension MxJT for models with multiple primary periods T.
}
\section{Methods}{
\describe{
\item{\code{signature(object = "unmarkedFit")}}{A fitted model object}
\item{\code{signature(object = "unmarkedFitDS")}}{A fitted model object}
\item{\code{signature(object = "unmarkedFitMPois")}}{A fitted model object}
\item{\code{signature(object = "unmarkedFitGMM")}}{A fitted model object}
\item{\code{signature(object = "unmarkedFitOccuCOP")}}{With \code{unmarkedFitOccuCOP} the object of a model fitted with \code{occuCOP}. Returns a matrix of \eqn{\lambda} the detection rate.}
}}
\keyword{methods}
