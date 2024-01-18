# unmarked 1.4.0

* Added count-data occupancy model (occuCOP)
* Added multi-scale occupancy model (goccu) 
* Added ZIP support to gdistsamp, gmultmix, and gpcount
* Fixed bug in TMB engine for occu that resulted in incorrect detection coefficient estimates when there were many interspersed NAs in the encounter history

# unmarked 1.3.3

* Increase required R version to 4.0
  
# unmarked 1.3.2

* Modernize some Cpp code to pass new LTO checks

# unmarked 1.3.1

* Remove log.grad function to pass CRAN checks

# unmarked 1.3.0

* Add support for terra package rasters
* Add plotEffects function for plotting marginal effects
* Better default names in fitLists
* Optional Shiny app for power analysis
* parboot now more robust to errors
* Add back temporarily removed occuMulti and colext vignettes
* Remove dependency on plyr package and move methods to imports
* Expand powerAnalysis vignette
* Many small bugfixes

# unmarked 1.2.4

* Convert vignettes to use rmarkdown
* Handle suggested packages in vignettes
* Remove occuMulti vignette due to AHMbook being temporarily off CRAN

# unmarked 1.2.3

* Add gdistremoval function to fit distance/removal models, see Amundson et al. 2014
* Add power analysis tools (powerAnalysis)
* Simulate datasets from scratch for any model type using simulate()
* Add penalized likelihood option to occuMulti, see Clipp et al. 2021
* Experimental random effects support for distsamp, multinomPois, and gdistremoval using TMB
* Improvements to predict() speed and better error messages
* Add vignettes for occuMulti, power analysis, simulation, and random effects
* Overhaul package tests and move to testthat
* New package website using pkgdown
* Move raster package from imports to suggests
* Fix assorted compilation warnings with newer versions of compilers on CRAN
* Remove call in TMB code to deprecated DOUBLE_XMIN
* Many bugfixes

# unmarked 1.1.1

* Fix address sanitizer problems with multmixOpen

# unmarked 1.1.0

* Add nmixTTD fitting function
* Add experimental random effects support and TMB engine to occu and pcount
* Add openMP support to some fitting functions (occuRN, gdistsamp, gmultmix, gpcount) for calculating likelihood in parallel
* Define STRICT_R_HEADERS in C++ code for compatibility with future Rcpp update
* Many bugfixes mainly related to predict()

# unmarked 1.0.1

* Fix LTO mismatches
* Automatically convert characters to factors in unmarkedFrames
* Many bugfixes, mainly related to predict()

# unmarked 1.0.0

* New functions 'distsampOpen' and 'multmixOpen' - open population versions of distsamp/gdistsamp and multinomPois/gmultmix
* Add 'predict' method for output from 'ranef', for generating posterior samples of the random effects and running a function on them
* Predict now correctly handles formulas containing functions and newdata with invalid factor levels
* Remove reshape2 dependency
* Bugfixes

# unmarked 0.13-1

* Fixes for compatibility with R 4.0

# unmarked 0.13-0

* New 'occuMS' function added for fitting multi-state occupancy models (single-season and dynamic)
* New 'occuTTD' function for fitting continuous time-to-detection occupancy models (single season and dynamic). Thanks to Jonathan Cohen for help with this
* New 'crossVal' function for doing cross-validation on fitted unmarked models and fitLists
* New 'vif' function for calculating variance inflation factors for fitted unmarked models
* Add ability to use complimentary-log-log link function in occu
* Add built-in dependent double observer pi function
* New C++ engines for gmultmix, gdistsamp, multinomPois, occuRN
* Approximate integrals in C++ engines with trapezoidal rule function instead of using Rdqags
* Misc minor bugfixes


# unmarked 0.12-3

* New 'occuMulti' function added by Ken Kellner

# unmarked 0.12-0

* Fixed mistake in turnover calculations in colext vignette (thanks to Giacomo Tavecchia) 
* added pcount.spHDS from AHM book. 
* updated predict method for pcount to include ZIP model
* Adam Smith added some parallel capabilities to the parboot functionality
* Adam Smith fixed formatMult conversion of julian date to factor
* Auriel Fournier fixed formatDistData to pad data with NA
* fixed error in obsToY for custom pi function

# unmarked 0.11-0

* Andy Royle is the new maintainer  
* Added Rebecca Hutchinson's penalized likelihood function occuPEN (experimental) 
* fixed bug in gmultmix to accommodate mixed sampling protocols (NA in count frequency vector is not counted in the constraint that multinomial cell probabilities sum to 1)  
* Changed variable 'ufp' to 'ufc' in ovenbird data and related functions.
* Removed constraint in pcountOpen that lambdaformula==omegaformula 
* Fixed bug in gdistsamp that caused error when NAs were present in half-normal model  
* Fixed bug in ranef (it was giving an error message for pcountOpen with the new dynamics options (Ricker and Gompertz) and working incorrectly for pcountOpen with immigration)
* Fixed bug in pcountOpen that occurred when covariates were time varying but not varying among sites 

# unmarked 0.10-6

* Fixed bug in C++ code that was causing problems on Solaris 

# unmarked 0.10-5

* Added new models of population dynamics to pcountOpen. Most changes contributed by Jeff Hostetler.  

# unmarked 0.10-4

* Added importFrom("plyr", "ldply") to NAMESPACE because "reshape" no longer depends on "plyr"

# unmarked 0.10-3

* RcppArmadillo was moved from "Depends" section of DESCRIPTION file to "LinkingTo"

# unmarked 0.10-2

* Thanks for Dirk Eddelbuettel for patch to deal with change in Armadillo's in-place reshape function. Serious problems might occur if you use a recent version of RcppArmadillo and an old version of unmarked.
* Dave Miller added another NA handling fix in occuFP(). I forgot to add this one in the previous version.

# unmarked 0.10-1

* Doc fixes requested by CRAN

# unmarked 0.10-0

* Fixed NA handling in occuFP()
* Fixed integration setting in C++ code that were causing segfaults when calling distsamp
* Replace raster:::layerNames() with raster:::names()
* distsamp() and gdistsamp() should be faster and more stable for some keyfun/survey combinations
