# Several vignettes (colext, power) take too long for CRAN to run
# Thus we have to pre-generate the results.
# The raw files (without results yet) are .Rmd.orig files
# These are ignored in package building
# The final files (with results) are .Rmd - these are the files that actually build the vignettes
# To generate an .Rmd from an .Rmd.orig (eg after updating relevant code)

knitr::knit("colext.Rmd.orig", output="colext.Rmd")
knitr::knit("powerAnalysis.Rmd.orig", output="powerAnalysis.Rmd")

# This will run all the R code in the .Rmd.orig file and save the results
# directly into the corresponding .Rmd file, which will then compile instantly on CRAN
# Note that this will also create some figure png files which should not be deleted
