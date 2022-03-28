context("gpcount fitting function")

test_that("unmarkedFrameGPC subset works",{
    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedFrameGPC(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)

    dat <- as(umf1, "data.frame")

    umf1.site1 <- umf1[1,]
    expect_equal(umf1.site1@y, y[1,, drop=FALSE])
    expect_equal(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    expect_equivalent(unlist(umf1.site1@obsCovs), oc$x3[1,])
    expect_equivalent(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    expect_equal(umf1.site1@numPrimary, 3)

    umf1.sites1and3 <- umf1[c(1,3),]

    expect_is(umf1.site1, "unmarkedFrameGPC")

    umf1.sites1and1 <- umf1[c(1,1),]

    umf1.obs1and2 <- umf1[,c(1,2)]

    expect_equivalent(dim(getY(umf1.obs1and2)), c(3,6))
    expect_equivalent(dim(siteCovs(umf1.obs1and2)), c(3,1))
    expect_equivalent(dim(obsCovs(umf1.obs1and2)), c(18,1))

    umf1.sites1and2.obs1and2 <- umf1[c(1,2),c(1,2)]
    expect_equal(class(umf1.sites1and2.obs1and2)[1], "unmarkedFrameGPC")
    expect_equivalent(dim(getY(umf1.sites1and2.obs1and2)), c(2,6))
    expect_equivalent(dim(siteCovs(umf1.sites1and2.obs1and2)), c(2,1))
    expect_equivalent(dim(obsCovs(umf1.sites1and2.obs1and2)), c(12,1))

    # THis doesn't work
    umf1.sites1and1.obs1and1 <- umf1[c(1,1),c(1,1)]
})

test_that("gpcount function works", {
  y <- matrix(c(0,0,0, 1,0,1, 2,2,2,
                3,2,3, 2,2,2, 1,1,1,
                NA,0,0, 0,0,0, 0,0,0,
                3,3,3, 3,1,3, 2,2,1,
                0,0,0, 0,0,0, 0,0,0), 5, 9, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,-1,4,-1))
  obsCovs <- list(o1 = matrix(seq(-3, 3, length=length(y)), 5, 9))
  obsCovs$o1[5,4:6] <- NA
  yrSiteCovs <- list(yr=matrix(c('1','2','2'), 5, 3, byrow=TRUE))
  yrSiteCovs$yr[4,2] <- NA

  expect_warning(umf <- unmarkedFrameGPC(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, numPrimary=3))

  expect_warning(fm <- gpcount(~x, ~yr, ~o1, data = umf, K=23))
  expect_equal(fm@sitesRemoved, integer(0))
  expect_equivalent(coef(fm),
        c(1.14754541, 0.44499137, -1.52079283, -0.08881542,
          2.52037155, -0.10950615), tol = 1e-5)

  # Check methods
  expect_warning(gp <- getP(fm))
  expect_equal(dim(gp), dim(y))

  expect_warning(pr <- predict(fm, 'lambda'))
  expect_equal(dim(pr), c(nrow(y), 4))

  nd <- data.frame(x=c(0,1))
  pr <- predict(fm, 'lambda', newdata=nd)
  expect_equal(dim(pr), c(2,4))
  expect_equal(pr[1,1], c(3.15045), tol=1e-4)

  res <- residuals(fm)
  expect_equal(dim(res), dim(y))

  expect_warning(r <- ranef(fm))
  expect_equal(dim(r@post), c(nrow(y), 24, 1))
  expect_equal(bup(r), c(7.31, 12.63, 1.30, 16.12, 2.04), tol=1e-3)

  expect_warning(s <- simulate(fm, 2))
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(y))

  expect_warning(pb <- parboot(fm, nsim=1))
  expect_is(pb, "parboot")

  # Check error when random effect in formula
  expect_error(gpcount(~(1|dummy),~1,~1,umf))

})

test_that("gpcount R and C++ engines give same results",{

  y <- matrix(c(0,0,0, 1,0,1, 2,2,2,
                3,2,3, 2,2,2, 1,1,1,
                NA,0,0, 0,0,0, 0,0,0,
                3,3,3, 3,1,3, 2,2,1,
                0,0,0, 0,0,0, 0,0,0), 5, 9, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,-1,4,-1))
  obsCovs <- list(o1 = matrix(seq(-3, 3, length=length(y)), 5, 9))
  yrSiteCovs <- list(yr=matrix(c('1','2','2'), 5, 3, byrow=TRUE))


  expect_warning(umf <- unmarkedFrameGPC(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, numPrimary=3))

  fm <- gpcount(~x, ~yr, ~o1, data = umf, K=23, control=list(maxit=1))
  fmR <- gpcount(~x, ~yr, ~o1, data = umf, K=23, engine="R", control=list(maxit=1))
  expect_equal(coef(fm), coef(fmR))
})
