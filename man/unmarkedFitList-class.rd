\name{unmarkedFitList-class}
\Rdversion{1.1}
\docType{class}
\alias{unmarkedFitList-class}
\alias{modSel,unmarkedFitList-method}
\alias{summary,unmarkedFitList-method}
\alias{coef,unmarkedFitList-method}
\alias{SE,unmarkedFitList-method}

\title{Class "unmarkedFitList"}
\description{Class to hold multiple fitted models from one of
	\code{unmarked}'s fitting functions}
\section{Objects from the Class}{
Objects can be created by using the \code{\link{fitList}} function.
}
\section{Slots}{
	 \describe{
    \item{\code{fits}:}{A \code{"list"} of models.}
  }
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "unmarkedFitList")}:
        Extract coefficients}
    \item{SE}{\code{signature(object = "unmarkedFitList")}:
        Extract standard errors}
    \item{modSel}{\code{signature(object = "unmarkedFitList")}:
		Model selection}
    \item{predict}{\code{signature(object = "unmarkedFitList")}:
		Model-averaged prediction}
	 }
}
\note{Model-averaging regression coefficients is intentionally not implemented.}
\seealso{
	\code{\link{fitList}},
	\code{\linkS4class{unmarkedFit}}
}
\examples{
showClass("unmarkedFitList")

data(linetran)
(dbreaksLine <- c(0, 5, 10, 15, 20))
lengths <- linetran$Length * 1000

ltUMF <- with(linetran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
	siteCovs = data.frame(Length, area, habitat), dist.breaks = dbreaksLine,
	tlength = lengths, survey = "line", unitsIn = "m")
	})

fm1 <- distsamp(~ 1 ~1, ltUMF)
fm2 <- distsamp(~ area ~1, ltUMF)
fm3 <- distsamp( ~ 1 ~area, ltUMF)

fl <- fitList(Null=fm1, A.=fm2, .A=fm3)
fl

coef(fl)
SE(fl)

ms <- modSel(fl, nullmod="Null")
ms


}
\keyword{classes}
