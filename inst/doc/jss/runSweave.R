library(cacheSweave)
Sweave("jss617.Rnw", driver=cacheSweaveDriver)
tools:::texi2dvi("jss617.tex", pdf=TRUE)

