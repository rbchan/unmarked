context("parboot")

y <- matrix(rep(0:1,10)[1:10],5,2)
siteCovs <- data.frame(x = c(0,2,3,4,1))
obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
fm <- occu(~ o1 + o2 ~ x, data = umf)

fitstats <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2,na.rm=TRUE)
    chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
    out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
    return(out)
}

test_that("parboot works", {
  pb <- parboot(fm, fitstats, nsim=3)
  expect_equal(dim(pb@t.star), c(3,3))

  # check show
  pb_out <- capture.output(pb)
  expect_equal(pb_out[4], "Parametric Bootstrap Statistics:")

  # check plot
  pdf(NULL)
  pl <- plot(pb)
  dev.off()
  expect_equal(pl, NULL)

  # check that report works
  rep_output <- capture.output(parboot(fm, fitstats, nsim=3, report=TRUE))
  expect_equal(substr(rep_output[1], 1,2), "t0")
})

test_that("parboot works in parallel",{
  skip_on_cran()
  skip_on_ci()
  # check parallel
  pb <- parboot(fm, nsim=101, parallel=TRUE, ncores=2)
  expect_equal(length(pb@t.star), 101)
})
