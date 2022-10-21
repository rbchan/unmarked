context("occuRN fitting function")
skip_on_cran()

test_that("occuRN can fit models",{

  data(birds)
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)

  # R and C engines give same result
  fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R", K=5, control=list(maxit=1))
  fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C", K=5, control=list(maxit=1))
  expect_equal(fm_R@AIC, fm_C@AIC)

  # survey occasion-specific detection probabilities
  fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C", K=10)
  #fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R")

  # check that output matches
  #expect_equivalent(coef(fm_C),coef(fm_R),tol=1e-5)

  # check output is correct
  expect_equivalent(coef(fm_C),
    c(0.7921122,-1.8328867,0.4268205,-0.1442194,0.4634105,0.7787513,
      0.8008794,1.0569827,0.8048578,0.8779660,0.9374874,0.7064848),tol=1e-3)

  # check methods
  gp <- getP(fm_C)
  expect_equal(dim(gp), dim(woodthrushUMF@y))

  pr <- predict(fm_C, 'state')
  expect_equal(dim(pr), c(50,4))
  expect_equal(pr[1,1], 2.204779, tol=1e-4)

  pr <- predict(fm_C, 'det')
  expect_equal(dim(pr), c(550,4))
  expect_equal(pr[1,1], 0.13806, tol=1e-4)

  res <- residuals(fm_C)
  expect_equal(dim(res), dim(woodthrushUMF@y))
  expect_equal(res[1,1], 0.73757, tol=1e-4)
  r <- ranef(fm_C, K=10)
  expect_equal(dim(r@post), c(50,11,1))

  s <- simulate(fm_C, 2)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(woodthrushUMF@y))

  pb <- parboot(fm_C, nsim=1)
  expect_is(pb, "parboot")

  # check error if random effect in formula
  expect_error(occuRN(~(1|dummy)~1, umf))
})

test_that("occuRN can handle NAs",{

  data(birds)
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)

  #Remove one observation
  woodthrushUMF@y[1,1] <- NA

  fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C", K=10)
  #fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R")

  # check that output matches
  #expect_equivalent(coef(fm_C),coef(fm_R),tol=1e-5)

  # check output is correct
  expect_equivalent(coef(fm_C),
    c(0.793042, -1.902789, 0.494098, -0.074573, 0.53074, 0.845903,
    0.867936, 1.123959, 0.871912, 0.944917, 1.004499, 0.773679), tol=1e-3)

  #Remove entire site
  woodthrush.bin_na <- woodthrush.bin
  woodthrush.bin_na[1,] <- NA
  woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin_na)

  expect_warning(fm_C <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="C", K=10))
  #fm_R <- occuRN(~ obsNum ~ 1, woodthrushUMF, engine="R")

  # check that site was removed
  expect_equivalent(fm_C@sitesRemoved,1)

  # check that output matches
  #expect_equivalent(coef(fm_C),coef(fm_R),tol=1e-5)

  # check output is correct
  expect_equivalent(coef(fm_C),
    c(0.783066, -1.920232, 0.448369, -0.009701, 0.490085, 0.814767,
    0.837669, 1.097903, 0.842467, 0.916831, 0.976707, 0.740672), tol=1e-3)
})

