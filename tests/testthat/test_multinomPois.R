context("multinomPois fitting function")
skip_on_cran()

test_that("unmarkedFrameMPois can be constructed",{
    y <- matrix(c(
        5, 3, 2,
        3, 3, 1,
        2, 0, 0,
        0, 0, 0,
        0, 0, 0), nrow=5, ncol=3, byrow=TRUE)

    sc <- data.frame(x1 = c(NA, 2, 3, 4, 3))
    oc <- list(x2 = matrix(c(
        1, 1, 1,
        3, NA, 1,
        0, 0, 1,
        NA, NA, NA,
        NA, 1, 0), nrow=5, ncol=3, byrow=TRUE))

    umf1 <- unmarkedFrameMPois(y = y, siteCovs = sc, obsCovs = oc,
        type="removal")
    expect_is(umf1, "unmarkedFrameMPois")

    o2y <- diag(ncol(y))
    o2y[upper.tri(o2y)] <- 1
    expect_equal(obsToY(umf1), o2y)

    expect_error(umf2 <- unmarkedFrameMPois(y=y, siteCovs=sc, obsCovs=oc, type="double"))
    umf2 <- unmarkedFrameMPois(y=y, siteCovs=sc,
                               obsCovs=lapply(oc, function(x) x[,1:2]),
                               type="double")
    expect_is(umf2, "unmarkedFrameMPois")

    umf3 <- unmarkedFrameMPois(y=y[,1:2], siteCovs=sc,
                               obsCovs=lapply(oc, function(x) x[,1:2]),
                               type="depDouble")
    expect_is(umf3, "unmarkedFrameMPois")

    expect_error(umf4 <- unmarkedFrameMPois(y=y, siteCovs=sc, obsCovs=oc, type="fake"))

    # error when depDouble and >2 samples
    expect_error(unmarkedFrameMPois(y=y, siteCovs=sc,
                               obsCovs=lapply(oc, function(x) x[,1:2]),
                               type="depDouble"))


})

test_that("multinomPois can fit a removal model",{

    y <- matrix(c(
        5, 3, 2,
        3, 3, 1,
        2, 0, 0,
        0, 0, 0,
        0, 0, 0), nrow=5, ncol=3, byrow=TRUE)

    sc <- data.frame(x1 = c(NA, 2, 3, 4, 3))
    oc <- list(x2 = matrix(c(
        1, 1, 1,
        3, NA, 1,
        0, 0, 1,
        NA, NA, NA,
        NA, 1, 0), nrow=5, ncol=3, byrow=TRUE))

    umf1 <- unmarkedFrameMPois(y = y, siteCovs = sc, obsCovs = oc,
        type="removal")

    #m1_R <- multinomPois(~1 ~1, umf1, engine="R")
    m1_C <- multinomPois(~1 ~1, umf1, engine="C")
    expect_equivalent(coef(m1_C), c(1.5257743, -0.2328092), tol=1e-5)
    #checkEqualsNumeric(coef(m1_R), coef(m1_C), tol=1e-5)

    #m2_R <- multinomPois(~x2 ~1, umf1, engine="R")
    expect_warning(m2_C <- multinomPois(~x2 ~1, umf1, engine="C"))
    expect_equivalent(coef(m2_C), c(1.9159845, 0.2248897, -0.1808144), tol=1e-5)
    expect_equal(m2_C@sitesRemoved, 4:5)
    #checkEqualsNumeric(coef(m2_R),coef(m2_C), tol=1e-5)

    #m3_R <- multinomPois(~x2 ~x1, umf1, engine="R")
    expect_warning(m3_C <- multinomPois(~x2 ~x1, umf1, engine="C"))
    expect_equivalent(m3_C@sitesRemoved, c(1, 4:5))
    expect_equivalent(coef(m3_C),
        c(1.9118525, -0.4071202, 8.3569943, 0.3232485), tol=1e-5)
    #checkEqualsNumeric(coef(m3_R),coef(m3_C), tol=1e-5)

    # check methods
    expect_warning(gp <- getP(m2_C))
    expect_equal(dim(gp), c(3,3))
    expect_equal(gp[1,1], 0.51101, tol=1e-4)

    expect_warning(pr <- predict(m2_C, 'state'))
    expect_equal(dim(pr), c(3,4))
    expect_equal(pr[1,1], 6.7936, tol=1e-4)

    nd <- data.frame(x2=c(0,1))
    pr <- predict(m2_C, 'det', newdata=nd)
    expect_equal(dim(pr), c(2,4))
    expect_equal(pr[1,1], 0.55598, tol=1e-4)

    res <- residuals(m2_C)
    expect_equal(dim(res), dim(umf1@y))

    expect_warning(r <- ranef(m2_C, K=50))
    expect_equal(dim(r@post), c(3,51,1))
    expect_equal(bup(r), c(10.794,6.9317,2.655), tol=1e-4)

    umf2 <- unmarkedFrameMPois(y=y, siteCovs=data.frame(x1=rnorm(5)), type="removal")
    m4 <- multinomPois(~1~x1, umf2)
    r <- ranef(m4, K=30)
    expect_equal(dim(r@post), c(5,31,1))

    expect_warning(s <- simulate(m2_C, 2, na.rm=FALSE))
    expect_equal(length(s), 2)

    expect_equal(dim(s[[1]]), dim(umf1@y))

    expect_warning(pb <- parboot(m2_C, nsim=1))
    expect_is(pb, "parboot")
})

test_that("multinomPois can fit a double observer model",{
  nSites <- 50
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  set.seed(9023)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }

  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)

  umf <- expect_warning(unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="double"))

  fm_C <- multinomPois(~observer-1 ~1, umf, engine="C")
  expect_equivalent(coef(fm_C), c(2.2586622, 0.1739752, -0.5685933), tol = 1e-5)
  expect_is(ranef(fm_C, K=30), "unmarkedRanef")

})

test_that("multinomPois can fit a dependent double observer model",{
  nSites <- 50
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1, p2*(1-p1))
  set.seed(9023)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 2)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:2]
  }
  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  umf <- expect_warning(unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="depDouble"))

  fm_C <- multinomPois(~observer-1 ~1, umf, engine="C")
  expect_equivalent(coef(fm_C), c(2.0416086, 0.7430343, 0.4564236), tol = 1e-5)
  expect_warning(r <- ranef(fm_C, K=30))
  expect_is(r, "unmarkedRanef")
})


test_that("multinomPois handles NAs",{
    y <- matrix(c(
        1, 0, 0,
        2, 1, 0,
        1, 0, 1,
        2, 1, 2,
        1, 0, 3,
        1, 1, 1), nrow=6, ncol=3, byrow=TRUE)
    oc <- matrix(c(
        1, 0,
        2, 1,
        1, 1,
        NA, 0,
        1, NA,
        NA, NA), nrow=6, ncol=2, byrow=TRUE)

    umf <- unmarkedFrameMPois(y = y, obsCovs = list(x=oc), type="double")

    expect_warning(m2 <- multinomPois(~x ~1, umf, starts=c(1.3, 0, 0.2)))
    expect_equal(m2@sitesRemoved, 4:6)

})

test_that("multinomPois can fit models with random effects",{
  set.seed(9023)
  nSites <- 50
  lambda <- 10
  p1 <- 0.5
  p2 <- 0.3
  cp <- c(p1*(1-p2), p2*(1-p1), p1*p2)
  N <- rpois(nSites, lambda)
  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }

  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  expect_warning(umf <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
            type="double"))
  fm <- multinomPois(~observer-1 ~1, umf)
  expect_true(inherits(fm, "unmarkedFitMPois"))
  expect_true(is.null(fm@TMB))
  pr <- predict(fm, "state")
  expect_equivalent(dim(pr), c(50,4))

  set.seed(1)
  nSites <- 100
  lambda <- 5
  sc <- data.frame(ref=sample(letters[1:10], nSites, replace=T),
                   x1=rnorm(nSites))
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)

  ef <- rnorm(10, 0, 0.4)
  names(ef) <- letters[1:10]
  lambda <- exp(log(lambda) + ef[sc$ref])
  N <- rpois(nSites, lambda)

  y <- matrix(NA, nSites, 3)
  for(i in 1:nSites) {
    y[i,] <- rmultinom(1, N[i], c(cp, 1-sum(cp)))[1:3]
  }
  expect_warning(umf2 <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
         type="double", siteCovs=sc))

  fm <- multinomPois(~observer-1 ~x1 + (1|ref), umf2)

  expect_true(inherits(fm@TMB, "list"))
  expect_equivalent(sigma(fm)$sigma, 0.3655, tol=1e-3)
  expect_true(inherits(randomTerms(fm), "data.frame"))
  pr <- predict(fm, type='state')
  pr2 <- predict(fm, "state", newdata=umf2@siteCovs[1:5,])
  expect_equivalent(dim(pr), c(100, 4))
  expect_equivalent(dim(pr2), c(5,4))

  # Make sure simulate accounts for random effects
  s <- simulate(fm, nsim=30)
  avg <- apply(sapply(s, function(x) x[,1]),1, mean)
  # average first count and predicted abundance should be highly correlated
  expect_true(cor(avg, pr$Predicted) > 0.7)

  umf2@y[1,1] <- NA
  umf2@y[2,] <- NA
  umf2@siteCovs$x1[3] <- NA
  umf2@obsCovs$observer[80] <- NA

  expect_warning(fm_na <- multinomPois(~observer-1 ~x1 + (1|ref), umf2))
  expect_true(inherits(fm_na, "unmarkedFitMPois"))

  expect_warning(umf3 <- unmarkedFrameMPois(y=y, obsCovs=list(observer=observer),
            piFun="fake", obsToY=umf@obsToY, siteCovs=sc))

  expect_error(multinomPois(~observer-1 ~x1 + (1|ref), umf3))

  # Site covs in detection formula
  expect_warning(fm <- multinomPois(~(1|ref)~1, umf2))
  expect_true(sigma(fm)$Model[1]=="p")
})
