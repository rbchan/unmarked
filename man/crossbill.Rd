\name{crossbill}
\alias{crossbill}
\docType{data}
\title{
  Detection/non-detection data on the European crossbill (\emph{Loxia
    curvirostra})
}
\description{
  267 1-kmsq quadrats were surveyed 3 times per year during 1999-2007.
}
\usage{data(crossbill)}
\format{
  A data frame with 267 observations on the following 58 variables.
  \describe{
    \item{\code{id}}{Plot ID}
    \item{\code{ele}}{Elevation}
    \item{\code{forest}}{Percent forest cover}
    \item{\code{surveys}}{a numeric vector}
    \item{\code{det991}}{Detection data for 1999, survey 1}
    \item{\code{det992}}{Detection data for 1999, survey 2}
    \item{\code{det993}}{Detection data for 1999, survey 3}
    \item{\code{det001}}{Detection data for 2000, survey 1}
    \item{\code{det002}}{a numeric vector}
    \item{\code{det003}}{a numeric vector}
    \item{\code{det011}}{a numeric vector}
    \item{\code{det012}}{a numeric vector}
    \item{\code{det013}}{a numeric vector}
    \item{\code{det021}}{a numeric vector}
    \item{\code{det022}}{a numeric vector}
    \item{\code{det023}}{a numeric vector}
    \item{\code{det031}}{a numeric vector}
    \item{\code{det032}}{a numeric vector}
    \item{\code{det033}}{a numeric vector}
    \item{\code{det041}}{a numeric vector}
    \item{\code{det042}}{a numeric vector}
    \item{\code{det043}}{a numeric vector}
    \item{\code{det051}}{a numeric vector}
    \item{\code{det052}}{a numeric vector}
    \item{\code{det053}}{a numeric vector}
    \item{\code{det061}}{a numeric vector}
    \item{\code{det062}}{a numeric vector}
    \item{\code{det063}}{Detection data for 2006, survey 3}
    \item{\code{det071}}{Detection data for 2007, survey 1}
    \item{\code{det072}}{Detection data for 2007, survey 2}
    \item{\code{det073}}{Detection data for 2007, survey 3}
    \item{\code{date991}}{Day of the season for 1999, survey 1}
    \item{\code{date992}}{Day of the season for 1999, survey 2}
    \item{\code{date993}}{Day of the season for 1999, survey 3}
    \item{\code{date001}}{Day of the season for 2000, survey 1}
    \item{\code{date002}}{a numeric vector}
    \item{\code{date003}}{a numeric vector}
    \item{\code{date011}}{a numeric vector}
    \item{\code{date012}}{a numeric vector}
    \item{\code{date013}}{a numeric vector}
    \item{\code{date021}}{a numeric vector}
    \item{\code{date022}}{a numeric vector}
    \item{\code{date023}}{a numeric vector}
    \item{\code{date031}}{a numeric vector}
    \item{\code{date032}}{a numeric vector}
    \item{\code{date033}}{a numeric vector}
    \item{\code{date041}}{a numeric vector}
    \item{\code{date042}}{a numeric vector}
    \item{\code{date043}}{a numeric vector}
    \item{\code{date051}}{a numeric vector}
    \item{\code{date052}}{a numeric vector}
    \item{\code{date053}}{a numeric vector}
    \item{\code{date061}}{a numeric vector}
    \item{\code{date062}}{a numeric vector}
    \item{\code{date063}}{a numeric vector}
    \item{\code{date071}}{a numeric vector}
    \item{\code{date072}}{a numeric vector}
    \item{\code{date073}}{Day of the season for 2007, survey 3}
  }
}
\source{
  Schmid, H. N. Zbinden, and V. Keller. 2004. Uberwachung der
    Bestandsentwicklung haufiger Brutvogel in der Schweiz,
    Swiss Ornithological Institute Sempach Switzerland
}
\seealso{\code{\link{Switzerland}} for corresponding covariate data
  defined for all 1-kmsq pixels in Switzerland. Useful for making
  species distribution maps.}
\examples{
data(crossbill)
str(crossbill)
}
\keyword{datasets}
