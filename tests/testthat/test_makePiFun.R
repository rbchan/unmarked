context("Generated piFuns")

test_that("getCombs returns all possible encounter histories",{
  # all-zero eh is removed

  out <- getCombs(1)
  expect_equal(out, 1)

  out <- getCombs(2)
  expect_equal(out, matrix(c(0,1,1,0,1,1), nrow=3, byrow=T))

  out <- getCombs(3)
  expect_equal(out, matrix(c(0,0,1,0,1,0,0,1,1,1,0,0,
                             1,0,1,1,1,0,1,1,1), nrow=7, byrow=T))

})

test_that("makeRemPiFun generates piFuns for varying-length removal periods",{

  out <- makeRemPiFun(c(0.5,1,1.5))
  expect_is(out, "function")
  p <- matrix(c(0.2,0.15,0.1,0.2,0.15,0.1), nrow=2, byrow=T)
  cp <- out(p)
  expect_is(cp, "matrix")
  expect_equal(dim(cp), c(2,3))
  expect_equal(cp[1], 0.1055728)

  p <- matrix(c(0.2,0.15,0.1,0.2), nrow=2, byrow=T)
  expect_error(out(p))

})

test_that("makeCrPiFun generates MRR pifun",{

  mrr <- makeCrPiFun(2)
  expect_is(mrr, "function")

  p <- matrix(c(0.2,0.15,0.1,0.2), nrow=2, byrow=T)

  cp <- mrr(p)
  expect_equal(cp, structure(c(0.12, 0.18, 0.17, 0.08, 0.03, 0.02),
                             .Dim = 2:3, .Dimnames = list(
                              NULL, c("01", "10", "11"))))
})

test_that("makeCrPiFunMb generates behavioral response pifun",{

  mrr <- makeCrPiFunMb(2)
  expect_is(mrr, "function")

  p <- matrix(c(0.2,0.15,0.1,0.2), nrow=2, byrow=T)
  cp <- mrr(p)
  expect_equal(cp, structure(c(0.16, 0.09, 0.17, 0.08, 0.03, 0.02),
                             .Dim = 2:3, .Dimnames = list(
                              NULL, c("01", "10", "11"))))

})

test_that("makeCrPiFunMh generates individ hetero pifun",{
  mrr <- makeCrPiFunMh(2)
  expect_is(mrr, "function")

  p <- matrix(c(0.2,0.15,0.1,0.2), nrow=2, byrow=T)
  cp <- mrr(p)
  expect_equal(cp, structure(c(0.16008585185222, 0.0906366734014539,
                               0.16008585185222,0.0906366734014539, 0.0413982898797996,
                               0.010483669158271), .Dim = 2:3, .Dimnames= list(
                                NULL, c("01", "10", "11"))))
})
