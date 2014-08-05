\name{mallard}
\alias{mallard}
\alias{mallard.y}
\alias{mallard.obs}
\alias{mallard.site}
\docType{data}
\title{Mallard count data}
\description{
Mallard repeated count data and covariates
}
\usage{data(mallard)}
\format{
A list with 3 components
	\describe{
	\item{mallard.y}{response matrix}
	\item{mallard.site}{site-specific covariates}
	\item{mallard.obs}{survey-specific covariates}
	}
}
\references{
Kery, M., Royle, J. A., and Schmid, H. (2005) Modeling Avaian Abundance from Replicated
 Counts Using Binomial Mixture Models. \emph{Ecological Applications} 15(4), 
 pp. 1450--1461.}
\examples{
data(mallard)
str(mallard.y)
str(mallard.site)
str(mallard.obs)
}
\keyword{datasets}
