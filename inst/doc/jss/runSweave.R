library(cacheSweave)
Sweave("unmarked.Rnw", driver=cacheSweaveDriver)
tools:::texi2dvi("unmarked.tex", pdf=TRUE)
