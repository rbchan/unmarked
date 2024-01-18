context("unmarkedFrames")

test_that("unmarkedFrame can be constructed",{
  M <- 10
  J <- 3
  y <- matrix(rbinom(J * M, 1, 0.5), M, J)
  siteCovs <- data.frame(a = rnorm(M), b = factor(gl(2,5)))
  umf <- unmarkedFrame(y = y, siteCovs = siteCovs)
  expect_is(umf, "unmarkedFrame")

  out <- capture.output(umf)
  expect_equal(out[1], "Data frame representation of unmarkedFrame object.")
  s <- capture.output(summary(umf))
  expect_equal(s[1], "unmarkedFrame Object")

  # convert to data frame
  df <- as(umf, "data.frame")
  expect_is(df, "data.frame")
})

test_that("obsToY works", {
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
    o2y <- obsToY(umf)

    expect_equal(o2y, matrix(1, 2, 3))
    oc.na <- is.na(oc)
    mult <- oc.na %*% o2y
    expect_is(mult, "matrix")
})

test_that("Characters are converted to factors when umf is constructed",{
  n <- 50   # number of sites
  T <- 4    # number of primary periods
  J <- 3    # number of secondary periods

  y <- matrix(0:1, n, J*T)

  #Site covs
  sc <- data.frame(x=rnorm(n), y=sample(letters, 50, replace=TRUE))
  expect_equal(sapply(sc, class), c(x="numeric", y="character"))
  expect_warning(umf <- unmarkedFrame(y, siteCovs=sc))

  umf <- expect_warning(unmarkedFrame(y, siteCovs=sc))
  expect_equal(sapply(siteCovs(umf), class), c(x="numeric", y="factor"))

  #Already factor
  sc2 <- data.frame(x=rnorm(n), y=factor(sample(letters, 50, replace=TRUE)))
  umf <- unmarkedFrame(y, siteCovs=sc2)
  expect_equal(sapply(siteCovs(umf), class), c(x="numeric", y="factor"))

  #Obs covs
  oc <- data.frame(x=rnorm(n*J*T), y=sample(letters, n*J*T, replace=TRUE))
  expect_equal(sapply(oc, class), c(x="numeric", y="character"))

  expect_warning(umf <- unmarkedFrame(y, obsCovs=oc))

  umf <- expect_warning(unmarkedFrame(y, obsCovs=oc))
  expect_equal(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  expect_true(is.null(siteCovs(umf)))

  #as list
  oc <- list(x=matrix(oc$x, nrow=n), y=matrix(oc$y, nrow=n))
  umf <- expect_warning(unmarkedFrameOccu(y, obsCovs=oc))
  expect_equal(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  expect_true(is.null(siteCovs(umf)))

  #Check conversion
  df <- as(umf, "data.frame")
  expect_equivalent(dim(df), c(50,36))

  #Yearly site covs
  ysc <- list(x=matrix(rnorm(n*T), nrow=n),
             y=matrix(sample(letters, n*T, replace=TRUE), nrow=n))
  umf <- expect_warning(unmarkedMultFrame(y, yearlySiteCovs=ysc, numPrimary=T))
  expect_equal(sapply(yearlySiteCovs(umf), class), c(x="numeric", y="factor"))
  expect_true(is.null(siteCovs(umf)))

  #All

  umf <- expect_warning(unmarkedMultFrame(y, yearlySiteCovs=ysc, obsCovs=oc,
                          siteCovs=sc, numPrimary=T))
  expect_equal(sapply(yearlySiteCovs(umf), class), c(x="numeric", y="factor"))
  expect_equal(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))
  expect_equal(sapply(obsCovs(umf), class), c(x="numeric", y="factor"))

  df <- as(umf, "data.frame")
  expect_equivalent(dim(df), c(50,46))
})

test_that("unmarkedMultFrame can be constructed",{
    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)
    expect_is(umf1, "unmarkedMultFrame")
    out <- capture.output(umf1)
    expect_equal(out[1], "Data frame representation of unmarkedFrame object.")

    s <- capture.output(summary(umf1))
    expect_equal(s[6], "Number of primary survey periods: 3 ")
})

test_that("unmarkedMultFrame subset works",{
    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)

    dat <- as(umf1, "data.frame")

    umf1.obs1 <- umf1[,1]
    expect_equal(umf1.obs1@y, y[,1:3])
    expect_equal(umf1.obs1@siteCovs, sc)
    expect_equivalent(unlist(umf1.obs1@obsCovs),
                       as.numeric(t(oc[[1]][,1:3])))
    expect_equivalent(unlist(umf1.obs1@yearlySiteCovs), ysc[[1]][,1])
    expect_equal(umf1.obs1@numPrimary, 1)

    umf1.obs1and3 <- umf1[,c(1,3)]

    umf1.site1 <- umf1[1,]
    expect_equal(umf1.site1@y, y[1,, drop=FALSE])
    expect_equal(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    expect_equivalent(unlist(umf1.site1@obsCovs), oc$x3[1,])
    expect_equivalent(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    expect_equal(umf1.site1@numPrimary, 3)

    umf1.sites1and3 <- umf1[c(1,3),]
})

test_that("unmmarkedMultFrame handles unequal secondary periods",{
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

  umf1 <- unmarkedMultFrame(y=y, numPrimary=4)
  expect_true(inherits(umf1, "unmarkedMultFrame"))

  expect_error(unmarkedMultFrame(y=y[,-1], numPrimary=4))
})

test_that("yearlySiteCovs processing works",{

  n <- 50   # number of sites
  T <- 4    # number of primary periods
  J <- 3    # number of secondary periods

  site <- 1:50
  years <- data.frame(matrix(rep(2010:2013, each=n), n, T))
  years <- data.frame(lapply(years, as.factor))
  dummy <- matrix(rep(c('a','b','c','d'),n),nrow=n,byrow=T)
  occasions <- data.frame(matrix(rep(1:(J*T), each=n), n, J*T))
  y <- matrix(0:1, n, J*T)

  umf <- expect_warning(unmarkedMultFrame(y=y,
        siteCovs = data.frame(site=site),
        obsCovs=list(occasion=occasions),
        yearlySiteCovs=list(year=years,dummy=dummy),
        numPrimary=T))

  as_df <- as(umf,'data.frame')

  expect_equivalent(dim(as_df),c(50,33))
  expect_true(all(names(as_df)[13:22] == c('site','year.1','year.2','year.3',
                                     'year.4','dummy.1','dummy.2','dummy.3',
                                     'dummy.4','occasion.1')))
  expect_true(all(as_df$year.1==2010))
  expect_true(all(as_df$dummy.1=='a'))


  umf2 <- unmarkedMultFrame(y=y,
        siteCovs = data.frame(site=site),
        obsCovs=list(occasion=occasions),
        numPrimary=T)

  as_df2 <- as(umf2,'data.frame')

  expect_equivalent(dim(as_df2),c(50,25))
})

test_that("lists provided to obsCovs or yearlySiteCovs must be named", {
    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)
    expect_is(umf1, "unmarkedMultFrame")


    oc <- list(matrix(1:27, 3))
    umf1 <- expect_error(unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3))

    oc <- list(x3 = matrix(1:27, 3))
    ysc <- list(matrix(1:9, 3))
    umf1 <- expect_error(unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3))

    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3), matrix(1:27, 3))

    umf1 <- expect_error(unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3))
})

test_that("covsToDF", {
  expect_equal(covsToDF(NULL, "obsCovs", 2, 3), NULL)
  
  df <- data.frame(x = rnorm(6), y = rnorm(6))
  expect_equal(covsToDF(df, "obsCovs", 2, 3),
               df)
  expect_error(covsToDF(df, "obsCovs", 2, 2))

  cl <- list(x = matrix(rnorm(6), 2, 3), y =matrix(rnorm(6), 2, 3))
  df_cl <- as.data.frame(lapply(cl, function(x) as.vector(t(x))))
  expect_equal(covsToDF(cl, "obsCovs", 3, 2),
               df_cl)
  expect_error(covsToDF(cl, "obsCovs", 2, 3))
})
