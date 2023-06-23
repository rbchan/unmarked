context("fitList and modSel methods")

skip_on_cran()

test_that("fitLists can be constructed",{
  y <- matrix(rep(1, 10), 5, 2)
  umf <- unmarkedFrameOccu(y = y, siteCovs=data.frame(x=-2:2),
        obsCovs= data.frame(z=-5:4))
  obsCovs(umf)[3, 1] <- NA
  fm1 <- occu(~ 1 ~ 1, data = umf)
  fm2 <- occu(~ 1 ~ x, data = umf)

  fits1.1 <- fitList(m1=fm1, m2=fm2)
  expect_equal(names(fits1.1@fits), c("m1","m2"))
  expect_message(fits1.2 <- fitList(fm1, fm2))
  expect_equal(names(fits1.2@fits), c("fm1","fm2"))
  fits2.1 <- fitList(fits = list(m1=fm1, m2=fm2))
  expect_equal(names(fits2.1@fits), c("m1","m2"))
  expect_message(fits2.2 <- fitList(fits = list(fm1, fm2)))
  expect_equal(names(fits2.2@fits), c("1","2"))

  expect_equal(fits1.1, fits2.1)

  expect_error(fitList(fm1, fm2, fits=list(fm1, fm2)))

  siteCovs(umf) <- data.frame(x=-3:1)
  fm2 <- occu(~ 1 ~ x, data = umf)
  expect_error(expect_warning(fitList(fm1, fm2)))   # Different umf used

  expect_warning(fm3 <- occu(~ z ~ 1, data = umf))
  expect_error(expect_warning(fitList(fm1, fm3)))   # Missing value problem
})

test_that("modSel method works",{
  y <- matrix(rep(1, 10), 5, 2)
  umf <- unmarkedFrameOccu(y = y, siteCovs=data.frame(x=-2:2),
        obsCovs= data.frame(z=-5:4))
  fm1 <- occu(~ 1 ~ 1, data = umf)
  fm2 <- occu(~ 1 ~ x, data = umf)

  fits <- fitList(m1=fm1, m2=fm2)
  ms1 <- modSel(fits)

  expect_true(all(is.na(ms1@Full$Rsq)))
  expect_equal(sum(ms1@Full$AICwt), 1)
  expect_equal(ms1@Full$delta[1L], 0)

  expect_error(modSel(fits, nullmod=fm2))

  ms2 <- modSel(fits, nullmod='m1')

  expect_equal(
        ms1@Full[,-which(colnames(ms1@Full)=="Rsq")],
        ms1@Full[,-which(colnames(ms2@Full)=="Rsq")]
  )

  # Fake hessian problem
  fm1@opt$hessian[] <- NA
  fm1@estimates@estimates$state@covMat[] <- NA
  fits2 <- fitList(m1=fm1, m2=fm2)
  ms3 <- modSel(fits2)
  expect_equal(coef(ms1), coef(ms3))

})

