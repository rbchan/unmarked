\name{modSel}
\title{Model selection results from an unmarkedFitList}
\description{Model selection results from an unmarkedFitList}
\alias{modSel}
\alias{modSel-methods}
\alias{unmarkedModSel-class}
\alias{show,unmarkedModSel-method}
\alias{summary,unmarkedModSel-method}
\alias{coef,unmarkedModSel-method}
\alias{SE,unmarkedModSel-method}
\arguments{
	\item{object}{an object of class "unmarkedFitList" created by the function
	  \code{\link{fitList}}.}
	\item{nullmod}{optional character naming which model in the
	\code{fitList} contains results from the null model. Only used in
		calculation of Nagelkerke's R-squared index.
	      }
}
\value{
A S4 object with the following slots
	\item{Full}{data.frame with formula, estimates, standard errors and model 
		selection information. Converge is optim convergence code. CondNum is
		model condition number. n is the number of sites. delta is delta AIC.
		cumltvWt is cumulative AIC weight.
        Rsq is Nagelkerke's (1991) R-squared index, which is only returned when 
        the nullmod argument is specified. }
	\item{Names}{matrix referencing column names of estimates (row 1) and 
	   standard errors (row 2).}
}
\note{Two requirements exist to conduct AIC-based model-selection and model-averaging in unmarked. First, the data objects (ie, unmarkedFrames) must be identical among fitted models. Second, the response matrix must be identical among fitted models after missing values have been removed. This means that if a response value was removed in one model due to missingness, it needs to be removed from all models.
} 
\references{Nagelkerke, N.J.D. (2004) A Note on a General Definition of the 
	Coefficient of Determination. \emph{Biometrika} 78, pp. 691-692.}
\author{Richard Chandler \email{rbchan@uga.edu}}
\examples{
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

ms <- modSel(fl, nullmod="Null")
ms

coef(ms)                            # Estimates only
SE(ms)                              # Standard errors only
(toExport <- as(ms, "data.frame"))  # Everything

}
