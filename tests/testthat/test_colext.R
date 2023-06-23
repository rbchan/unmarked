context("colext fitting function")
skip_on_cran()

# Simulate data
set.seed(123)
nsites <- 6
nyr <- 4
nrep <- 2
y <- matrix(c(
        1,0, 1,1, 0,0, 0,0,
        1,1, 0,0, 0,0, 0,0,
        0,0, 0,0, 0,0, 0,0,
        0,0, 1,1, 0,0, 0,0,
        1,1, 1,0, 0,1, 0,0,
        0,0, 0,0, 0,0, 1,1), nrow=nsites, ncol=nyr*nrep, byrow=TRUE)

sc <- data.frame(sc1 = rnorm(nsites))
oc <- matrix(rnorm(nsites*nyr*nrep), nsites, nyr*nrep)
ysc <- matrix(rnorm(nsites*nyr), nsites, nyr)

test_that("unmarkedMultFrame construction works",{

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, yearlySiteCovs=list(ysc=ysc),
                            obsCovs=list(oc=oc), numPrimary=4)
  expect_is(umf1, "unmarkedMultFrame")

  expect_error(unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc[1:5,]), numPrimary=4))

  expect_error(unmarkedMultFrame(y=y, siteCovs=sc, yearlySiteCovs=list(ysc=ysc),
                            obsCovs=list(oc=oc), numPrimary=3))

  plot(umf1)
})


test_that("colext model fitting works", {

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)

  fm1 <- colext(~1, ~1, ~1, ~1, umf1)
    expect_equivalent(coef(fm1),
        c(0.1422577, -1.4950576,  0.2100365,  1.1998444),
        tol=1e-6)

  # With site covs
  fm2 <- colext(~sc1, ~1, ~1, ~1, umf1)
  expect_equivalent(coef(fm2), c(1.3423, -6.2788,-1.5831,0.1413,1.1638),
                    tol=1e-4)

  # With obs covs
  fm3 <- colext(~1, ~1, ~1, ~oc, umf1)
  expect_equivalent(coef(fm3),
        c(0.1433,-1.4975,0.2082,1.2002,-0.03822),
        tol=1e-4)

  # With yearly site covs
  fm4 <- colext(~1, ~ysc, ~ysc, ~1, umf1)
  expect_equivalent(coef(fm4),
                    c(0.2662,-2.0534,-1.0579,0.2165,0.6877,1.10342), tol=1e-4)

  # ranef
  r <- ranef(fm4)
  expect_equal(dim(r@post), c(nsites, nrep, nyr))
  expect_equal(dim(bup(r)), c(nsites, nyr))

})

test_that("colext handles missing values",{

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)

  umf2 <- umf1
  umf2@y[1,3] <- NA

  fm1 <- colext(~1, ~1, ~1, ~1, umf2)
  expect_is(fm1, "unmarkedFitColExt")

  umf3 <- umf1
  umf3@y[1,] <- NA
  expect_warning(fm2 <- colext(~1, ~1, ~1, ~1, umf3))
  expect_is(fm2, "unmarkedFitColExt")
  expect_equal(fm2@sitesRemoved, 1)

  umf4 <- umf1
  umf4@y[1,3:4] <- NA
  fm3 <- colext(~1, ~1, ~1, ~1, umf4)
  expect_is(fm3, "unmarkedFitColExt")

  umf5 <- umf1
  umf5@siteCovs$sc1[1] <- NA
  expect_warning(fm4 <- colext(~sc1, ~1, ~1, ~1, umf5))
  expect_warning(pr <- predict(fm4, 'psi'))
  expect_equal(nrow(pr), 5)

umf5 <- umf1
  umf5@obsCovs$oc[1] <- NA
  expect_warning(fm4 <- colext(~1, ~1, ~1, ~oc, umf5))
  expect_warning(pr <- predict(fm4, 'det'))
  expect_equal(nrow(pr), nsites*nyr*nrep)
  expect_true(all(is.na(pr[1,])))

  umf5 <- umf1
  umf5@yearlySiteCovs$ysc[1] <- NA
  # This should work, right?
  expect_error(expect_warning(fm4 <- colext(~1, ~1, ~ysc, ~1, umf5)))

})

test_that("colext errors when random effects are in formula",{
  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)
  expect_error(colext(~(1|dummy), ~1, ~ysc, ~1, umf1))
})

test_that("colext methods work",{

  umf1 <- unmarkedMultFrame(y=y, siteCovs=sc, obsCovs=list(oc=oc),
                            yearlySiteCovs=list(ysc=ysc), numPrimary=nyr)
  fm1 <- colext(~sc1, ~1, ~ysc, ~oc, umf1)

  pdf(NULL)
  plot(fm1)
  dev.off()
  res <- residuals(fm1)
  expect_equal(dim(res), c(6,8))
  r <- ranef(fm1)
  expect_equal(dim(r@post), c(nsites, nrep, nyr))
  pr1 <- predict(fm1, 'psi')
  expect_equal(nrow(pr1), 6)
  pr2 <- predict(fm1, 'col')
  expect_equal(nrow(pr2), nsites*nyr)
  pr3 <- predict(fm1, 'det')
  expect_equal(nrow(pr3), nsites*nyr*nrep)

  nd <- data.frame(sc1=c(0,1))
  pr4 <- predict(fm1, 'psi', newdata=nd)
  expect_equal(nrow(pr4), 2)
})
