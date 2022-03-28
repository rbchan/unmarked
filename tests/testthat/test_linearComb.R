context("linearComb and backTransform")

y <- matrix(rep(0:1,10)[1:10],5,2)
siteCovs <- data.frame(x = c(0,2,3,4,1))
obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
fm <- occu(~ o1 + o2 ~ x, data = umf)

lc <- linearComb(fm, type='state', c(1,0.1))

test_that("linearComb works",{

  expect_is(lc, "unmarkedLinComb")
  expect_equal(lc@estimate, as.numeric(c(1,0.1) %*% coef(fm, 'state')))
  out <- capture.output(lc)
  expect_equal(out[1], "Linear combination(s) of Occupancy estimate(s)")

  df <- as(lc, "data.frame")
  expect_is(df, "data.frame")
  expect_equal(df[1,1], lc@estimate)
})

test_that("backTransform works",{

  bt <- backTransform(lc)
  expect_is(bt, "unmarkedBackTrans")
  out <- capture.output(bt)
  expect_equal(out[1], "Backtransformed linear combination(s) of Occupancy estimate(s)")

  df <- as(bt, "data.frame")
  expect_is(df, "data.frame")
  expect_equal(df[1,1], 0.9998549)
})
