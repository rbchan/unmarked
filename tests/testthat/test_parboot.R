context("parboot")
skip_on_cran()

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

fitstats2 <- function(fm, na.rm=TRUE) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    resids <- residuals(fm)
    sse <- sum(resids^2, na.rm=na.rm)
    chisq <- sum((observed - expected)^2 / expected, na.rm=na.rm)
    freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=na.rm)
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

  # check that report argument gives warning
  expect_warning(parboot(fm, fitstats, nsim=3, report=TRUE))
})

test_that("parboot works in parallel",{
  skip_on_cran()
  skip_on_ci()
  # check parallel
  pb <- parboot(fm, nsim=10, parallel=TRUE, ncores=2)
  expect_equal(length(pb@t.star), 10)
})

test_that("parboot handles failing model fits", {

  fail_func <- function(x){
    rand <- rnorm(1)
    if(rand > 0.5){
      stop("fail")
    }
    return(rand)
  }

  set.seed(123)
  expect_warning(pb <- parboot(fm, nsim=20, statistic=fail_func))
  expect_equal(nrow(pb@t.star), 13)

  # Error message when all parboot samples are bad

  # force error only when running function on new simulated datasets,
  # but not for original dataset
  fail_func <- function(x){
    if(round(x@AIC, 5) == 23.29768){
      return(0)
    } else {
      stop("fail")
    }
  }

  set.seed(123)
  expect_error(pb2 <- parboot(fm, nsim=20, statistic=fail_func))
})

test_that("parboot handles failing model fits in parallel", {
  skip_on_cran()
  skip_on_ci()
  fail_func <- function(x){
    rand <- rnorm(1)
    if(rand > 0.5){
      stop("fail")
    }
    return(rand)
  }

  set.seed(123)
  expect_warning(pb <- parboot(fm, nsim=20, statistic=fail_func, parallel=TRUE))
  expect_true(nrow(pb@t.star) < 20)

})

test_that("parboot handles statistic functions with additional arguments", {

  opt_func <- function(x, y){
    res <- mean(residuals(x), na.rm=TRUE)
    c(res=res, y=y)
  }

  pb <- parboot(fm, nsim=10, statistic=opt_func, y=0.1)
  expect_equal(colnames(pb@t.star), c("res", "y"))
  expect_true(all(pb@t.star[,"y"]==0.1))
})

test_that("parboot handles statistic functions with additional arguments in parallel", {
  skip_on_cran()
  skip_on_ci()
  opt_func <- function(x, y){
    res <- mean(residuals(x), na.rm=TRUE)
    c(res=res, y=y)
  }

  pb <- parboot(fm, nsim=10, statistic=opt_func, y=0.1, parallel=TRUE)
  expect_equal(colnames(pb@t.star), c("res", "y"))
  expect_true(all(pb@t.star[,"y"]==0.1))
})

# To show the fitstats function in parboot example works when there are NAs
test_that("fitstats2 function handles NAs", {
  umf2 <- umf
  umf2@y[1,1] <- NA
  fm2 <- occu(~ o1 + o2 ~ x, data = umf2)
  expect_error(expect_warning(show(parboot(fm2, nsim=3, statistic=fitstats2, na.rm=FALSE))))
  pb <- parboot(fm2, nsim=3, statistic=fitstats2, na.rm=TRUE)
  expect_is(pb, "parboot")
  expect_equal(nrow(pb@t.star), 3)
})
