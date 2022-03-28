context("pcountOpen fitting function")

test_that("unmarkedFramePCO subset works",{
    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedFramePCO(
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
    expect_is(umf1.site1, "unmarkedFramePCO")

    umf1.sites1and3 <- umf1[c(1,3),]

    # subset by primary period
    umf_sub <- umf1[,1:2]
    expect_equal(umf_sub@numPrimary, 2)
})

test_that("pcountOpen can fit a null model",{
  y <- matrix(c(
      3, 2, 1, 4,
      3, 4, 2, 1,
      0, 1, 2, 3,
      5, 3, 3, 4,
      2, 4, 3, 3), 5, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:20)
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10,
      starts=c(1, 0, 0, 7))
  expect_equivalent(coef(fm1),
                     c(0.9565311, 0.2741022, 0.1352888, 7.0041290),
                     tol = 1e-5)

  fm2 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="gamma",
                    K=10)
  expect_equivalent(coef(fm2), c(1.8219354, 8.7416638, -0.2873611),
      tol = 1e-3)

  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="omega",
                    K=30)
  expect_equivalent(coef(fm3), c(1.8861091, -1.3102890, -0.4934883),
      tol = 1e-5)

  # methods
  gp <- getP(fm1)
  expect_equal(dim(gp), dim(umf@y))

  res <- residuals(fm1)
  expect_equal(dim(res), dim(umf@y))
  expect_equal(res[1,1], 0.39970, tol=1e-4)

  s <- simulate(fm1, 2)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(umf@y))

  r <- ranef(fm1)
  expect_equal(dim(r@post), c(5,11,4))
  expect_equal(dim(bup(r)), c(5,4))
  expect_equal(bup(r)[1], 3.00236, tol=1e-5)

  pb <- parboot(fm1, nsim=1)
  expect_is(pb, "parboot")

  # Check error when random effect in formula
  expect_error(pcountOpen(~(1|dummy),~1,~1,~1, data=umf))

})


test_that("pcountOpen handles NAs",{

  y1 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA,
      2, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf1 <- unmarkedFramePCO(y = y1, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf1, se=FALSE, K=10,
      starts=c(1.6, 0.24, 1.16, -0.268))
  expect_equivalent(coef(fm1),
      c(0.9536852, 0.4611643, -0.8655834, 3.2154420), tol = 1e-5)

  y2 <- matrix(c(
      1, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1,
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  oc <- y1 + -2:3

  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- list(o1 = oc)
  ysc <- list(o2=oc)
  umf2 <- unmarkedFramePCO(y = y2, siteCovs = siteCovs, obsCovs = obsCovs,
      yearlySiteCovs=ysc, numPrimary=4)

  expect_warning(fm2.1 <- pcountOpen(~1, ~1, ~1, ~o1, data = umf2, se=FALSE, K=10,
      starts=c(1.4, -1.3, 1.8, -1.1, 0.7)))
  expect_equivalent(coef(fm2.1),
      c(1.239636, -2.085200, 1.770919, -0.602612, 1.255386),
                     tol = 1e-4)

  expect_warning(fm2.2 <- pcountOpen(~1, ~1, ~o2, ~1, data = umf2, se=FALSE, K=10,
      starts=c(1.2, -1, 2, 0, 0)))
  expect_equivalent(coef(fm2.2),
      c(1.3242059, 0.8439311, -2.8217070, -10.1414153, 0.1176959),
      tol = 1e-5)

  expect_warning(fm2.3 <- pcountOpen(~1, ~o2, ~1, ~1, data = umf2, se=FALSE, K=10,
      starts=c(1, 0, 0, -5, -1)))
  expect_equivalent(coef(fm2.3),
      c(0.7013386, 0.5277811, -0.2350951, -1.8346326, 4.7771974),
                     tol = 1e-2)

  y3 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA,
      NA, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf3 <- unmarkedFramePCO(y = y3, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  expect_warning(fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf3, se=FALSE, K=10,
      starts=c(1.5, 0, 1, 0)))
  expect_equivalent(coef(fm3),
                     c(0.9714456, 0.4481042, -0.8920462, 4.0379739 ),
                     tol = 1e-4)
  expect_equal(fm3@sitesRemoved, 6)



  y4 <- matrix(c(
      NA, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1,
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  go4 <- matrix(c(
      NA, NA, NA, # remove y[1, 2:4]
      NA, 1, 2,   # remove y[2, 2]
      0, NA, 2,   # remove y[3, 3]. Creates an interior NA
      5, 1, NA,   # remove y[4, 4]. Creates an end NA
      1, NA, NA,  # remove y[5, 3:4]
      NA, NA, 1), 6, 3, byrow=TRUE)
  o2y <- matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1), 3, 3, byrow=TRUE)
  y4.na <- is.na(go4) %*% o2y
  y4.2 <- y4
  y4.2[,-1][y4.na>0] <- NA
  #y4.2

  umf4 <- unmarkedFramePCO(y=y4, yearlySiteCovs=list(go4=cbind(go4, 1)),
      numPrimary=4)

  expect_warning(fm4.1 <- pcountOpen(~1, ~go4, ~1, ~1, umf4, se=FALSE,
      starts=c(.8, .5, -.3, -1.5, 6)))
  expect_equal(fm4.1@sitesRemoved, 1)

  expect_warning(fm4.2 <- pcountOpen(~1, ~1, ~go4, ~1, umf4, se=FALSE,
      starts=c(.8, 0, 5, -5, 7)))
  expect_equal(fm4.2@sitesRemoved, 1)



    # Now with secondary sampling periods

    y5 <- matrix(c(
        2,2,  2,2,  1,0,
        3,2,  1,1,  2,2,
        0,0,  1,1,  1,1,
        3,3,  2,0,  3,3), 4, 6, byrow=TRUE)

    umf5 <- unmarkedFramePCO(y=y5, numPrimary=3)
    fm5 <- pcountOpen(~1, ~1, ~1, ~1, umf5, se=FALSE, K=10)
    expect_equivalent(coef(fm5),
        c(0.7269958, -0.3484145, 0.1494188, 1.9391898), tol=1e-5)

    y6 <- y5
    y6[1,1] <- y6[2,3:4] <- y6[3,5:6] <- y6[4,6] <- NA
    umf6 <- unmarkedFramePCO(y=y6, numPrimary=3)
    fm6 <- pcountOpen(~1, ~1, ~1, ~1, umf6, se=FALSE, K=10)
    expect_equivalent(coef(fm6),
        c(0.7945817, -0.4340502, 0.5614526, 1.4161393), tol=1e-5)

    y7 <- y5
    oc7 <- y6 + -2:1
    umf7 <- unmarkedFramePCO(y=y7, obsCovs=list(oc=oc7), numPrimary=3)
    expect_warning(fm7 <- pcountOpen(~1, ~1, ~1, ~oc, umf7, se=FALSE, K=10))
    expect_equivalent(coef(fm7),
        c(1.1986029, -9.9298367, 12.5760064, -0.8876606, 0.9525870),
                       tol=1e-4)

    y8 <- y5
    ysc8 <- matrix(1:3, #rnorm(12),
                   4, 3, byrow=TRUE)
    ysc8[2,1] <- NA # note this will fail if it's ysc8[1,1] <- NA
    umf8 <- unmarkedFramePCO(y=y8, yearlySiteCovs=list(ysc=ysc8),
                             numPrimary=3)
  ## !!!!!!!!!!!
    expect_warning(fm8 <- pcountOpen(~1, ~1, ~ysc, ~1, umf8, se=FALSE, K=10))
    expect_equivalent(coef(fm8),
        c(0.7278796, -0.8770411, 0.9170578, 0.0399341, 1.8956210),
                       tol=1e-4)

})


test_that("pcountOpen handles different lengths of time between primary periods (delta)",{
    M <- 5
    T <- 4
    y <- matrix(c(
        NA, 2, 1, 4,
        3, NA, 2, 1,
        0, 1, NA, 3,
        5, 3, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    if(!exists("formatDelta"))
        formatDelta <- unmarked:::formatDelta
    dates <- matrix(c(1,3,5,7), M, T, byrow=TRUE)
    delta <- formatDelta(dates, is.na(y))
    ans <- matrix(c(
        1, 2, 2, 2,
        1, 2, 4, 2,
        1, 2, 2, 4,
        1, 2, 2, 2,
        1, 2, 2, 2), M, T, byrow=TRUE)

    expect_equal(delta, ans)

    dates2 <- matrix(c(
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8), M, T, byrow=TRUE)
    delta2 <- formatDelta(dates2, is.na(y))
    ans2 <- matrix(c(
        2, 3, 2, 2,
        1, 3, 5, 2,
        2, 2, 2, 4,
        1, 3, 2, 2,
        2, 3, 2, 2), M, T, byrow=TRUE)

    expect_equal(delta2, ans2)

    dates3 <- matrix(as.integer(c(
      2, NA, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8)), M, T, byrow=TRUE)
    expect_error(unmarkedFramePCO(y=y, primaryPeriod=dates3,
                                    numPrimary=4))

    dates4 <- dates2
#    dates4[is.na(y)] <- NA
    mode(dates4) <- "integer"
    delta4 <- formatDelta(dates4, is.na(y))
    umf <- unmarkedFramePCO(y=y, primaryPeriod=dates4, numPrimary=4)
    fm <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10, starts=c(1.2, 0, 1.4, 1.2))
    expect_equivalent(coef(fm),
        c(1.35779948, 0.11911809, -0.06946651, 5.78090618),
#        c(1.2543989, -0.5429887, 0.6715887, 5.6593500),
    tol = 1e-5)

    y5 <- matrix(c(
        1, NA, 1, 4,
        NA, 3, 2, 1,
        0, 1, 2, NA,
        NA, NA, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    dates5 <- matrix(c(
        2, NA, 6, 8,
        NA, 4, 6, 8,
        2, 4, 6, NA,
        NA, NA, 6, NA,
        2, 4, 6, 8), M, T, byrow=TRUE)
    ans5 <- matrix(c(
        1, NA, 4, 2,
        NA, 2, 2, 2,
        1, 2, 2, NA,
        NA, NA, 4, NA, # 4 not 5 b/c primary period 1 is day 2
        1, 2, 2, 2), M, T, byrow=TRUE)
    delta5 <- formatDelta(dates5, is.na(y5))
    expect_equal(delta5, ans5)

    dates6 <- y6 <- matrix(c(2L, 1L), 1, 2)
    expect_error(unmarkedFramePCO(y=y6, primaryPeriod=dates6,
                                    numPrimary=2))

})

test_that("pcountOpen can fit models with various dynamics",{
    set.seed(3)
    M <- 20
    T <- 5
    lambda <- 4
    gamma <- .2
    omega <- 0.7
    p <- 0.7
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega)
	G[,t] <- rpois(M, gamma*N[,t])
	N[,t+1] <- S[,t] + G[,t]
    }
    y[] <- rbinom(M*T, N, p)
    #colMeans(y)
    umf <- unmarkedFramePCO(y = y, numPrimary=T)

    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="autoreg")
    expect_equivalent(coef(m1),
                       c(1.5457081, -1.2129776,  0.5668830,  0.4987492),
                       tolerance=1e-5)

    m2 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="notrend")
    expect_equivalent(coef(m2),
                       c(1.2131713,  0.7301736,  1.1949289),
                       tolerance=1e-5)

    m3 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, dynamics="trend")
    expect_equivalent(coef(m3),
                       c(1.67211946, -0.06534021, 0.18287762),
                       tolerance=1e-5)

})


test_that("pcountOpen simulate and fitted methods work",{

     set.seed(3)
     M <- 50
     T <- 5
     lambda <- 4
     gamma <- 1.5
     omega <- 0.8
     p <- 0.7
     y <- N <- matrix(NA, M, T)
     S <- G <- matrix(NA, M, T-1)
     N[,1] <- rpois(M, lambda)
     for(t in 1:(T-1)) {
             S[,t] <- rbinom(M, N[,t], omega)
             G[,t] <- rpois(M, gamma)
             N[,t+1] <- S[,t] + G[,t]
             }
     y[] <- rbinom(M*T, N, p)

     umf <- unmarkedFramePCO(y = y, numPrimary=T)

     m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20)
     m1_sim <- simulate(m1)
     expect_equivalent(m1_sim[[1]][1,],c(3,4,4,2,2))

     m2 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, fix='gamma')
     m2_sim <- simulate(m2)
     expect_equivalent(m2_sim[[1]][1,],c(4,4,6,5,3))
     m2_fit <- fitted(m2)
     expect_equivalent(m2_fit[1,1],3.448029,tol=1e-4)

     m3 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20, fix='omega')
     m3_sim <- simulate(m3)
     expect_equivalent(m3_sim[[1]][1,],c(2,1,3,5,4))
     m3_fit <- fitted(m3)
     expect_equivalent(m3_fit[1,1], 2.481839, tol=1e-4)
})


test_that("pcountOpen predict method works",{
     set.seed(3)
     M <- 10
     T <- 3
     lambda <- 4
     gamma <- 1.5
     omega <- 0.8
     p <- 0.7
     y <- N <- matrix(NA, M, T)
     S <- G <- matrix(NA, M, T-1)
     N[,1] <- rpois(M, lambda)
     for(t in 1:(T-1)) {
             S[,t] <- rbinom(M, N[,t], omega)
             G[,t] <- rpois(M, gamma)
             N[,t+1] <- S[,t] + G[,t]
             }
     y[] <- rbinom(M*T, N, p)

     sc <- data.frame(x1=rnorm(M), x2 = runif(M, 0, 1),
                 x3=rnorm(M), x4=rnorm(M))

     umf <- unmarkedFramePCO(y = y, numPrimary=T, siteCovs=sc)

     m1 <- pcountOpen(~x1, ~x2, ~x3, ~x4, umf, K=9)

     #Make sure predicting with newdata works
     p1 <- predict(m1, type = "lambda",
            newdata = data.frame(x1 = rnorm(5), x2 = rnorm(5)))
     p2 <- predict(m1, type = "gamma",
            newdata = data.frame(x1 = rnorm(5), x2 = rnorm(5)))
     p3 <- predict(m1, type = "omega",
            newdata = data.frame(x3 = rnorm(5), x4 = rnorm(5)))
     p4 <- predict(m1, type = "det",
            newdata = data.frame(x3 = rnorm(5), x4 = rnorm(5)))

     are_df <- sapply(list(p1,p2,p3,p4), function(x) inherits(x, "data.frame"))
     expect_true(all(are_df))

     m2 <- pcountOpen(~x1, ~x2, ~x3, ~x4, umf, K=9, mixture="ZIP")
     expect_warning(pr <- predict(m2,'lambda'))
     expect_is(pr, "data.frame")
     expect_equivalent(as.numeric(pr[1,]), c(2.572621,NA,NA,NA), tol=1e-4)
})

