context("distsamp fitting function")


y <- matrix(rep(4:1, 10)[1:10], 5, 2, byrow=TRUE)
siteCovs <- data.frame(x = c(0, 2, 3, 4, 1))

test_that("unmarkedFrameDS identifies problems with inputs",{
  #Check error thrown when length(tlength!=nrow(y))
 expect_error(unmarkedFrameDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(0, 5, 10)/1000, survey="line", tlength=rep(1, (5-1)),
        unitsIn="km"))
  #Check error thrown when length(dist.breaks) != J+1
  expect_error(unmarkedFrameDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(5,10)/1000, survey="line", tlength=rep(1, 5),
        unitsIn="km"))

  #Check error when obs covs are provided
  oc <- data.frame(z=rnorm(10))
  expect_error(unmarkedFrameDS(y=y, siteCovs=siteCovs, obsCovs=oc))

  umf <- unmarkedFrameDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(0, 5, 10)/1000, survey="line", tlength=rep(1, 5),
        unitsIn="km")
  expect_is(umf, "unmarkedFrameDS")
  s <- capture.output(summary(umf))
  expect_equal(s[1], "unmarkedFrameDS Object")

  expect_error(obsCovs(umf) <- oc)

  # histogram
  pdf(NULL)
  pl <- hist(umf)
  expect_is(pl, "histogram")
  dev.off()

})

test_that("distsamp works with covariates", {
  umf <- unmarkedFrameDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(0, 5, 10)/1000, survey="line", tlength=rep(1, 5),
        unitsIn="km")
  fm <- distsamp(~ x ~ x, data = umf)

  lam <- fm['state']
  det <- fm['det']

  expect_equivalent(coef(lam), c(1.4340999, -0.1102387), tolerance = 1e-4)
  expect_equivalent(coef(det), c(-4.64686395, -0.09337832), tolerance = 1e-4)

  lam.lc <- linearComb(fm, type = 'state', c(1, 2))
  det.lc <- linearComb(fm, type = 'det', c(1, 2))

  expect_equivalent(coef(lam.lc), 1.213623, tol = 1e-4)
  expect_equivalent(coef(det.lc), -4.833621, tol = 1e-4)

  expect_equivalent(coef(backTransform(lam.lc)), 3.365655, tol = 1e-4)
  expect_equivalent(coef(backTransform(det.lc)), 0.007957658, tol = 1e-4)

})

test_that("distsamp methods work",{

  umf <- unmarkedFrameDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(0, 5, 10)/1000, survey="line", tlength=rep(1, 5),
        unitsIn="km")
  fm <- distsamp(~ x ~ x, data = umf)

  pr <- predict(fm, 'state')
  expect_equal(dim(pr), c(5,4))
  expect_equal(pr[1,1], 4.19586, tol=1e-4)

  pr <- predict(fm, 'det')
  expect_equal(dim(pr), c(5,4))
  expect_equal(pr[1,1], 0.00959, tol=1e-4)

  nd <- data.frame(x=c(0,1))
  pr <- predict(fm, 'state', newdata=nd)
  expect_equal(dim(pr), c(2,4))

  res <- residuals(fm)
  expect_equal(dim(res), dim(y))
  expect_equal(res[1,1], -0.01333, tol=1e-4)

  r <- ranef(fm, K=50)
  expect_is(r, "unmarkedRanef")
  expect_equal(dim(r@post), c(5,51,1))

  expect_error(hist(fm))

  pdf(NULL)
  plot(fm)
  dev.off()
})

test_that("distsamp ranef method works",{

   set.seed(344)
    lambda <- 10
    sigma <- 20
    npts <- 10
    radius <- 50
    breaks <- seq(0, 50, by=10)
    A <- (2*radius)^2 / 10000 # Area (ha) of square containing circle
    y <- matrix(0, npts, length(breaks)-1)
    N <- integer(npts)
    for(i in 1:npts) {
        M <- rpois(1, lambda * A) # Individuals within the square
        xy <- cbind(x=runif(M, -radius, radius),
                    y=runif(M, -radius, radius))
        d <- apply(xy, 1, function(x) sqrt(x[1]^2 + x[2]^2))
        d <- d[d <= radius]
        N[i] <- length(d)
        if(length(d)) {
            p <- exp(-d^2 / (2 * sigma^2)) # half-normal
            d <- d[rbinom(length(d), 1, p) == 1]
            y[i,] <- table(cut(d, breaks, include.lowest=TRUE))
        }
    }

    umf1 <- unmarkedFrameDS(y = y, survey="point",
                            dist.breaks=breaks, unitsIn="m")
    m1 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)))
    m2 <- distsamp(~1 ~1, umf1, starts=c(log(5), log(20)),
                    output="abund")

    re1 <- ranef(m1, K=20)
    re2 <- ranef(m2, K=20)

    expect_equal(mode1 <- bup(re1, stat="mode"), bup(re2, "mode"))
    expect_equal(confint(re1), confint(re2))

    ar1 <- as(re1, "array")

  expect_equivalent(colSums(ar1), c(
    0.000000e+00, 2.334960e-01, 8.517322e-01, 1.524261e+00, 1.811577e+00,
    1.691348e+00, 1.421738e+00, 1.085003e+00, 7.119743e-01, 3.898376e-01,
    1.782052e-01, 6.895313e-02, 2.296231e-02, 6.685198e-03, 1.725009e-03,
    3.991224e-04, 8.362689e-05, 1.600128e-05, 2.816112e-06, 4.586885e-07,
    6.951721e-08), tolerance=1e-6)
})

test_that("distsamp line keyfunctions work",{
    y <- structure(c(7, 7, 12, 9, 9, 11, 9, 5, 7, 6, 25, 26, 30, 26, 23,
        24, 20, 33, 26, 32, 5, 3, 8, 7, 1, 4, 4, 7, 7, 6, 3, 1, 1, 4,
        4, 4, 3, 6, 2, 3), .Dim = c(10L, 4L))
    umf <- unmarkedFrameDS(y = y, dist.breaks=c(0, 3, 15, 18, 20),
        survey="line", unitsIn="m", tlength=rep(100, nrow(y)))

    fm.halfnorm <- distsamp(~1~1, umf)
    D <- backTransform(fm.halfnorm, type="state")
    S <- backTransform(fm.halfnorm, type="det")
    expect_equivalent(coef(D), 129.5509, tol=1e-4)
    expect_equivalent(SE(D), 9.446125, tol=1e-4)
    expect_equivalent(coef(S), 18.15386, tol=1e-4)
    expect_equivalent(SE(S), 2.893362, tol=1e-4)

    fm.exp <- distsamp(~1~1, umf, keyfun="exp", starts=c(4, 0))
    D <- backTransform(fm.exp, type="state")
    S <- backTransform(fm.exp, type="det")
    expect_equivalent(coef(D), 144.8802, tol=1e-4)
    expect_equivalent(SE(D), 14.31655, tol=1e-4)
    expect_equivalent(coef(S), 31.75738, tol=1e-4)
    expect_equivalent(SE(S), 9.711254, tol=1e-4)

    fm.haz <- distsamp(~1~1, umf, keyfun="hazard", starts=c(4, 3, 1))
    D <- backTransform(fm.haz, type="state")
    Sh <- backTransform(fm.haz, type="det")
    Sc <- backTransform(fm.haz, type="scale")
    expect_equivalent(coef(D), 137.0375, tol=1e-4)
    expect_equivalent(SE(D), 16.82505, tol=1e-4)
    expect_equivalent(coef(Sh), 15.90262, tol=1e-4)
    expect_equivalent(SE(Sh), 5.099981, tol=1e-4)
    expect_equivalent(coef(Sc), 0.8315524, tol=1e-4)
    expect_equivalent(SE(Sc), 0.4753275, tol=1e-4)

    fm.unif <- distsamp(~1~1, umf, keyfun="uniform")
    D <- backTransform(fm.unif, type="state")
    expect_equivalent(coef(D), 107.5000, tol=1e-4)

    expect_equivalent(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    expect_equivalent(coef(fm.exp),
                       coef(update(fm.exp, engine="R")))
    expect_equivalent(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    expect_equivalent(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    expect_equivalent(coef(fm.unif), coef(update(fm.unif, engine="R")))
    expect_equivalent(coef(fm.haz), coef(update(fm.haz, engine="R")), tol=1e-5)

})

test_that("distsamp point keyfunctions work",{
    y <- structure(c(1, 0, 0, 0, 0, 0, 3, 1, 1, 0, 16, 15, 18, 14, 22,
        24, 12, 20, 20, 21, 10, 9, 9, 5, 6, 6, 6, 9, 5, 6, 6, 6, 4, 2,
        6, 3, 3, 3, 1, 4), .Dim = c(10L, 4L))

    umf <- unmarkedFrameDS(y = y, dist.breaks=c(0, 3, 15, 18, 20),
        survey="point", unitsIn="m", tlength=rep(100, 20))

    fm.halfnorm <- distsamp(~1~1, umf)
    D <- backTransform(fm.halfnorm, type="state")
    S <- backTransform(fm.halfnorm, type="det")
    expect_equivalent(coef(D), 316.1711, tol=1e-4)
    expect_equivalent(SE(D), 37.08797, tol=1e-4)
    expect_equivalent(coef(S), 18.05958, tol=1e-4)
    expect_equivalent(SE(S), 3.341798, tol=1e-4)

    fm.exp <- distsamp(~1~1, umf, keyfun="exp", starts=c(6, 0))
    D <- backTransform(fm.exp, type="state")
    S <- backTransform(fm.exp, type="det")
    expect_equivalent(coef(D), 369.7526, tol=1e-4)
    expect_equivalent(SE(D), 68.11901, tol=1e-4)
    expect_equivalent(coef(S), 28.90848, tol=1e-4)
    expect_equivalent(SE(S), 11.66219, tol=1e-4)

    fm.haz <- distsamp(~1~1, umf, keyfun="hazard", starts=c(5, 3, 1))
    D <- backTransform(fm.haz, type="state")
    Sh <- backTransform(fm.haz, type="det")
    Sc <- backTransform(fm.haz, type="scale")
    expect_equivalent(coef(D), 266.3911, tol=1e-4)
    expect_equivalent(SE(D), 20.45144, tol=1e-4)
    expect_equivalent(coef(Sh), 18.69351, tol=1e-4)
    expect_equivalent(SE(Sh), 0.8950444, tol=1e-4)
    expect_equivalent(coef(Sc), 5.797366, tol=1e-4)
    expect_equivalent(SE(Sc), 4.054381, tol=1e-4)

    fm.unif <- distsamp(~1~1, umf, keyfun="uniform")
    D <- backTransform(fm.unif, type="state")
    expect_equivalent(coef(D), 236.3451, tol=1e-4)

    expect_equivalent(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    expect_equivalent(coef(fm.exp),
                       coef(update(fm.exp, engine="R")),tol=1e-5)
    expect_equivalent(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))
    expect_equivalent(coef(fm.halfnorm),
                       coef(update(fm.halfnorm, engine="R")))

})

test_that("getP works with distsamp",{

  data(issj)
  jayumf <- unmarkedFrameDS(y=as.matrix(
  issj[,1:3]),
  siteCovs=data.frame(scale(issj[,c("elevation","forest","chaparral")])),
  dist.breaks=c(0,100,200,300), unitsIn="m", survey="point")

  hn <- distsamp(~1 ~1, jayumf)
  neg <- distsamp(~1 ~1, jayumf,keyfun="exp")
  unif <- distsamp(~1 ~1, jayumf, keyfun="unif")
  haz <- distsamp(~1 ~1, jayumf[1:100,], keyfun="hazard")

  expect_equivalent(getP(hn)[1,], c(0.08634098, 0.09873522, 0.02369782),
                     tol=1e-5)
  expect_equivalent(getP(neg)[1,], c(0.1111111, 0.3333333, 0.5555556),
                     tol=1e-5)
  expect_equivalent(getP(unif)[1,], c(0.1111111, 0.3333333, 0.5555556),
                     tol=1e-5)
  expect_equivalent(getP(haz)[1,], c(0.0702,0.2107,0.35117),
                     tol=1e-3)
})

test_that("distsamp works with random effects",{

  data(linetran)
  umf <- unmarkedFrameDS(y=as.matrix(linetran[,1:4]), siteCovs=linetran[,6:7],
                         survey="line", tlength=linetran$Length, unitsIn='m',
                         dist.breaks=c(0,10,20,30,40))

  hn <- distsamp(~1~area+(1|habitat), umf)
  ex <- distsamp(~1~area+(1|habitat), umf, keyfun="exp")
  hz <- distsamp(~1~area+(1|habitat), umf, keyfun="hazard")
  un <- distsamp(~1~area+(1|habitat), umf, keyfun="uniform")
  mods <- list(hn=hn, ex=ex, hz=hz, un=un)
  expect_true(all(sapply(mods, function(x) is.list(x@TMB))))

  sigs <- sapply(mods, function(x) sigma(x)$sigma)
  expect_true(all(sigs < 0.01) & all(sigs > 0.0001))

  pr <- lapply(mods,  function(x) predict(x, "state"))
  expect_true(all(sapply(pr, inherits, "data.frame")))

  data(pointtran)
  umf <- unmarkedFrameDS(y=as.matrix(pointtran[,1:4]), siteCovs=pointtran[,6:7],
                         survey="point", unitsIn='m',
                         dist.breaks=c(0,10,20,30,40))

  hn <- distsamp(~1~area+(1|habitat), umf)
  ex <- distsamp(~1~area+(1|habitat), umf, keyfun="exp")
  hz <- distsamp(~1~area+(1|habitat), umf, keyfun="hazard")
  un <- distsamp(~1~area+(1|habitat), umf, keyfun="uniform")
  mods <- list(hn=hn, ex=ex, hz=hz, un=un)
  expect_true(all(sapply(mods, function(x) is.list(x@TMB))))

  sigs <- sapply(mods, function(x) sigma(x)$sigma)
  expect_true(all(sigs < 0.01) & all(sigs > 0.0001))

  pr <- lapply(mods,  function(x) predict(x, "state"))
  expect_true(all(sapply(pr, inherits, "data.frame")))

})
