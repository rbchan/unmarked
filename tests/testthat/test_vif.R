context("vif method")
skip_on_cran()

test_that("vif works with occu models",{

  set.seed(123)
  data(frogs)
  pferUMF <- unmarkedFrameOccu(pfer.bin)
  siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))
  #Add some correlated covariates
  obsvar2 = rnorm(numSites(pferUMF) * obsNum(pferUMF))
  obsvar3 = rnorm(numSites(pferUMF) * obsNum(pferUMF),mean=obsvar2,sd=0.5)
  obsCovs(pferUMF) <- data.frame(
                        obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)),
                        obsvar2=obsvar2,obsvar3=obsvar3)

  fm <- occu(~ obsvar1+obsvar2+obsvar3 ~ 1, pferUMF)

  #No type provided
  expect_error(vif_vals <- vif(fm))

  #Wrong type provided
  expect_error(vif_vals <- vif(fm, type='fake'))

  #Not enough covs
  expect_error(vif_vals <- vif(fm, type='state'))

  #Get values for det
  vif_vals <- vif(fm, type='det')
  expect_equivalent(vif_vals, c(1.002240,4.49552,4.490039),tol=1e-4)
  expect_equal(names(vif_vals),c('obsvar1','obsvar2','obsvar3'))

  #Compare to typical way of calculating
  set.seed(123)
  vt <- lm(obsvar2~obsvar1+obsvar3,data=obsCovs(pferUMF))
  vt_2 <- 1/(1-summary(vt)$r.squared)
  expect_equivalent(vif_vals[2],vt_2,tol=0.1)

})

test_that("vif works with multinomPois models",{

  set.seed(123)
  data(ovendata)
  ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
  siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])),
  type = "removal")
  fm1 <- multinomPois(~ 1 ~ ufc + trba, ovenFrame)

  vif_vals <- vif(fm1, type='state')

  expect_equivalent(vif_vals, c(1.285886,1.285886), tol=1e-4)

})
