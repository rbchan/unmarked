context("occuPEN fitting function")
skip_on_cran()

test_that("occuPEN can fit simple models",{

  y <- matrix(rep(1,10)[1:10],5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occuPEN(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  expect_equivalent(occ,1)

  det <- coef(backTransform(det))
  expect_equivalent(det,1)

  bt <- backTransform(fm, type = 'state')
  expect_equivalent(coef(bt), 1)

  bt <- backTransform(fm, type = 'det')
  expect_equivalent(coef(bt), 1)

  #Check error when random effect in formula
  expect_error(occuPEN(~(1|dummy)~1, umf))

  y <- matrix(rep(0,10)[1:10],5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occuPEN(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  expect_equivalent(occ, 0, tolerance = 1e-4)

  det <- coef(backTransform(det))
  expect_equivalent(det,0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'state')
  expect_equivalent(coef(bt), 0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'det')
  expect_equivalent(coef(bt), 0, tolerance = 1e-4)

})

test_that("occuPEN can fit models with covariates",{

  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occuPEN(~ o1 + o2 ~ x, data = umf)
  fmR <- occuPEN(~o1 +o2~x, data=umf, engine="R")
  expect_equal(coef(fm), coef(fmR))
  fm1 <- occuPEN(~ o1 + o2 ~ x, data = umf,lambda=1,pen.type="Bayes")
  fm2 <- occuPEN(~ o1 + o2 ~ x, data = umf,lambda=1,pen.type="Ridge")
  MPLEla <- computeMPLElambda(~ o1 + o2 ~ x, data = umf)
  fm3 <- occuPEN(~ o1 + o2 ~ x, data = umf,lambda=MPLEla,pen.type="MPLE")

  occ <- fm['state']
  det <- fm['det']

  occ1 <- fm1['state']
  det1 <- fm1['det']

  occ2 <- fm2['state']
  det2 <- fm2['det']

  occ3 <- fm3['state']
  det3 <- fm3['det']

  expect_error(occ <- coef(backTransform(occ)))

  expect_equivalent(coef(occ), c(8.590737, 2.472220), tolerance = 1e-4)
  expect_equivalent(coef(det), c(0.44457, -0.14706, 0.44103), tolerance = 1e-4)

  expect_equivalent(coef(occ1), c(0.7171743, 0.6977753), tolerance = 1e-4)
  expect_equivalent(coef(det1), c(0.08143832, -0.06451574,  0.28695210), tolerance = 1e-4)

  expect_equivalent(coef(occ2), c(1.009337e+01,  4.329662e-04), tolerance = 1e-4)
  expect_equivalent(coef(det2), c(0.25892308, -0.09459618,  0.31092107), tolerance = 1e-4)

  expect_equivalent(coef(occ3), c(8.590738, 2.472220), tolerance = 1e-4)
  expect_equivalent(coef(det3), c(0.4445733, -0.1470601,  0.4410251), tolerance = 1e-4)

  occ.lc <- linearComb(fm, type = 'state', c(1, 0.5))
  det.lc <- linearComb(fm, type = 'det', c(1, 0.3, -0.3))

  expect_equivalent(coef(occ.lc), 9.826848, tol = 1e-4)
  expect_equivalent(coef(det.lc), 0.2681477, tol = 1e-4)

  expect_equivalent(coef(backTransform(occ.lc)), 1, tol = 1e-4)
  expect_equivalent(coef(backTransform(det.lc)), 0.5666381, tol = 1e-4)

  expect_error(backTransform(fm, type = "state"))
  expect_error(backTransform(fm, type = "det"))

  fitted <- fitted(fm)
  expect_equivalent(fitted, structure(c(0.5738, 0.5014, 0.4318, 0.38581, 0.50171, 0.53764,
0.46563, 0.40283, 0.39986, 0.79928), .Dim = c(5L, 2L)), tol = 1e-5)

  expect_error(occuPEN_CV(~ o1 + o2 ~ x, data = umf, k=15))
  fmCV <- occuPEN_CV(~ o1 + o2 ~ x, data = umf)
  expect_equivalent(fmCV@chosenLambda, 1, tol = 1e-4)
  expect_equivalent(fmCV@lambdaScores, c(31.423777, 15.603297, 12.330360, 10.130768,  8.981720,  8.572523,  8.572841, 8.798436,  9.153270,  9.543802), tol = 1e-4)


  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  expect_equivalent(computeMPLElambda(~ o1 + o2 ~ x, data = umf),4.017164e-11)
  expect_error(fm <- occuPEN(~ o1 + o2 ~ x, data = umf,pen.type="none"))
  expect_error(fm <- occuPEN(~ o1 + o2 ~ 1, data = umf,pen.type="MPLE"))
  expect_error(fm <- occuPEN(~ 1 ~ 1, data = umf,pen.type="Ridge"))
  expect_error(fm <- occuPEN_CV(~ o1 + o2 ~ x, data = umf,lambda=c(0)))
  expect_error(fm <- occuPEN_CV(~ o1 + o2 ~ x, data = umf,foldAssignments=c(1,2,3,4,5),k=6))

  # nonparboot
  nbp <- nonparboot(fm, B=2)
  expect_is(nbp@covMatBS, "matrix")
  nbp_cv <- nonparboot(fmCV, B=2)
  expect_is(nbp_cv@covMatBS, "matrix")

})

test_that("occuPEN can handle NAs",{

  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  expect_warning(fm <- occuPEN(~ o1 + o2 ~ x, data = umf))
  expect_equal(fm@sitesRemoved, 3)
  expect_equivalent(coef(fm), c(8.70123, 4.58255, 0.66243, -0.22862, 0.58192), tol = 1e-5)

  obsCovs[10,2] <- NA
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  expect_warning(fm <- occuPEN(~ o1 + o2 ~ x, data = umf))
  expect_equal(fm@sitesRemoved, 3)
  expect_equivalent(coef(fm), c(8.91289, 1.89291, -1.42471, 0.67011, -8.44608), tol = 1e-5)

})


test_that("occuPEN can handle offsets",{

  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occuPEN(~ o1 + o2 ~ offset(x), data = umf)
  expect_equivalent(coef(fm),
                     structure(c(9.74361, 0.44327, -0.14683, 0.44085), .Names = c("psi(Int)",
"p(Int)", "p(o1)", "p(o2)")), tol = 1e-5)
  fm <- occuPEN(~ o1 + offset(o2) ~ offset(x), data = umf)
  expect_equivalent(coef(fm), structure(c(8.59459, 0.97574, -0.3096), .Names = c("psi(Int)",
"p(Int)", "p(o1)")), tol=1e-5)

})
