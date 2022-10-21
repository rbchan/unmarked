context("crossVal method")
skip_on_cran()

set.seed(123)
data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

fm <- occu(~ obsvar1 ~ 1, pferUMF[1:20,])


test_that("crossVal works with occu models",{
  set.seed(123)
  kfold <- crossVal(fm, method='Kfold', folds=10)

  expect_equal(nrow(kfold@stats),10)

  expect_equal(as.numeric(kfold@stats[1,]),
                c(0.3790110,0.3014053), tolerance=1e-4)

  holdout <- crossVal(fm, method='holdout', holdoutPct=0.25)

  expect_equal(as.numeric(holdout@stats[1,]),
                       c(0.296829,0.262929), tolerance=1e-4)

  leave <- crossVal(fm, method='leaveOneOut')

  expect_equal(nrow(leave@stats),20)
  expect_equal(as.numeric(leave@stats[1,]),
               c(0.5985,0.5012), tolerance=1e-4)

  show_output <- capture.output(leave)
  expect_equal(show_output[1], "Method: leave-one-out")

})

test_that("crossVal works in parallel",{
  skip_on_cran()
  skip_on_ci()

  set.seed(123)
  kfold <- crossVal(fm, method='Kfold', folds=10)
  set.seed(123)
  kfold_par <- crossVal(fm, method='Kfold', folds=10, parallel=TRUE, ncores=2)
  expect_equal(kfold@stats, kfold_par@stats)


})

test_that("custom statistics functions work",{

  expect_error(crossVal(fm, statistic=function(x) "fake"))

  new_stat <- function(object){
    c(mean_res = mean(residuals(object),na.rm=T))
  }

  kfold_custom <- crossVal(fm, statistic=new_stat)
  expect_equal(length(kfold_custom@stats[,1]), 10)
})

test_that("crossValList can be constructed",{

  fm <- occu(~ obsvar1 ~ 1, pferUMF[1:20,])
  fm2 <- occu(~1 ~1, pferUMF[1:20,])

  fl <- fitList(fm2=fm2,fm=fm)
  cvlist <- crossVal(fl, method='Kfold')

  expect_is(cvlist, "unmarkedCrossValList")
  expect_equal(length(cvlist@stats_list),2)
  show_output <- capture.output(cvlist)
  expect_equal(show_output[1], "Method: k-fold (10 folds)")
})

test_that("crossVal works with multinomPois",{

  set.seed(123)
  data(ovendata)
  ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
    type = "removal")
  fm1 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame[1:20,])

  mout <- crossVal(fm1, method='Kfold')
  expect_equal(as.numeric(mout@stats[1,]),
                       c(0.25859,0.17974), tolerance=1e-4)

})
