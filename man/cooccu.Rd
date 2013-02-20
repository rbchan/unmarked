\name{cooccu}
\alias{cooccu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Co-occurrence model
}
\description{
A co-occurrence model allowing for (possibly) dependent occurrence,
false negative errors, and false positive errors.
}
\usage{
cooccu(psiformulaA = ~1, psiformulaB = ~1, gammaformula = NULL, pformulaA = ~1, pformulaB = ~1, fpformulaA = NULL, fpformulaB = NULL, data, interaction = c("dominance", "general"), fpmodel = c("full", "confusion"), starts, method = "BFGS", se = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psiformulaA}{
%%     ~~Describe \code{psiformulaA} here~~
}
  \item{psiformulaB}{
%%     ~~Describe \code{psiformulaB} here~~
}
  \item{gammaformula}{
%%     ~~Describe \code{gammaformula} here~~
}
  \item{pformulaA}{
%%     ~~Describe \code{pformulaA} here~~
}
  \item{pformulaB}{
%%     ~~Describe \code{pformulaB} here~~
}
  \item{fpformulaA}{
%%     ~~Describe \code{fpformulaA} here~~
}
  \item{fpformulaB}{
%%     ~~Describe \code{fpformulaB} here~~
}
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{interaction}{
%%     ~~Describe \code{interaction} here~~
}
  \item{fpmodel}{
%%     ~~Describe \code{fpmodel} here~~
}
  \item{starts}{
%%     ~~Describe \code{starts} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{se}{
%%     ~~Describe \code{se} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Richard Chandler
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nSites <- 100
nReps <- 5
yA <- matrix(rbinom(nSites*nReps, 1, 0.1), nSites, nReps)
yB <- matrix(rbinom(nSites*nReps, 1, 0.1), nSites, nReps)
umf <- unmarkedFrameCo(yA=yA, yB=yB)
fm1 <- cooccu(psiformulaA=~1, psiformulaB=~1, #gammaformula=~1,
              data=umf)
fm1
fm2 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=~1,
              data=umf)
fm2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
