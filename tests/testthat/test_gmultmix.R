context("gmultmix fitting function")
skip_on_cran()

test_that("unmarkedFrameGMM construction works",{

  y <- matrix(0:3, 5, 4)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
  yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  expect_is(umf, "unmarkedFrameGMM")

  expect_error(unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="fake", numPrimary=2))

  # error when >2 sampling occasions per primary period and type depDouble
  y2 <- cbind(y, y[,1:2])
  expect_error(unmarkedFrameGMM(y = y2, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="depDouble", numPrimary=2))

})

test_that("unmarkedFrameGMM subset works",{
  y <- matrix(1:27, 3)
  sc <- data.frame(x1 = 1:3)
  ysc <- list(x2 = matrix(1:9, 3))
  oc <- list(x3 = matrix(1:27, 3))

  umf1 <- unmarkedFrameGMM(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3,
        type="removal")

  dat <- as(umf1, "data.frame")

  umf1.site1 <- umf1[1,]
  expect_equal(umf1.site1@y, y[1,, drop=FALSE])
  expect_equal(umf1.site1@siteCovs, sc[1,, drop=FALSE])
  expect_equivalent(unlist(umf1.site1@obsCovs), oc$x3[1,])
  expect_equivalent(unlist(umf1.site1@yearlySiteCovs),
                    ysc$x2[1,, drop=FALSE])
  expect_equal(umf1.site1@numPrimary, 3)

  umf1.sites1and3 <- umf1[c(1,3),]
  expect_equal(class(umf1.site1)[1], "unmarkedFrameGMM")

  umf1.sites1and1 <- umf1[c(1,1),]
  umf1.obs1and2 <- umf1[,c(1,2)]

  expect_equivalent(dim(getY(umf1.obs1and2)), c(3,6))
  expect_equivalent(dim(siteCovs(umf1.obs1and2)), c(3,1))
  expect_equivalent(dim(obsCovs(umf1.obs1and2)), c(18,1))

  umf1.sites1and2.obs1and2 <- umf1[c(1,2),c(1,2)]
  expect_equivalent(dim(getY(umf1.sites1and2.obs1and2)), c(2,6))
  expect_equivalent(dim(siteCovs(umf1.sites1and2.obs1and2)), c(2,1))
  expect_equivalent(dim(obsCovs(umf1.sites1and2.obs1and2)), c(12,1))

  # THis doesn't work
  umf1.sites1and1.obs1and1 <- umf1[c(1,1),c(1,1)]
})

test_that("gmultmix removal model works",{
  y <- matrix(0:3, 5, 4)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
  yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C"))

  expect_equal(fm_C@sitesRemoved, 3)
  coef_truth <- c(2.50638554, 0.06226627, 0.21787839, 6.46029769, -1.51885928,
            -0.03409375, 0.43424295)
  #checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
  expect_equivalent(coef(fm_C), coef_truth, tol = 1e-5)

  # NAs in obsCovs
  obsCovs[10,2] <- NA
  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C"))

  expect_equal(fm_C@sitesRemoved, 3)
  #checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
  expect_equivalent(coef(fm_C), coef_truth, tol = 1e-5)

  # NAs in ysc
  yrSiteCovs[2, 1] <- NA
  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23, engine="C"))

  coef_truth <- c(1.17280104, 0.37694710, 2.38249795, 2.87354955, -0.83875134,
            -0.08446507, 1.88056826)
  #checkEqualsNumeric(coef(fm_R), coef_truth, tol = 1e-5)
  expect_equivalent(coef(fm_C), coef_truth, tol = 1e-5)

  #Negative binomial
  #fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, mixture="NB", K=23, engine="R")
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, mixture="NB",
                                  K=23, engine="C"))
  expect_equivalent(coef(fm_C), c(1.1819, 0.3738,2.4571,4.3633,-0.8734,-0.08211,
                                  1.86049,9.38619), tol=1e-4)

  #Check methods
  expect_warning(gp <- getP(fm_C))
  expect_equal(dim(gp), c(4,4)) # missing site dropped

  expect_warning(pr <- predict(fm_C, 'lambda'))
  expect_equal(dim(pr), c(4,4))

  nd <- data.frame(x=c(0,1))
  pr <- predict(fm_C, 'lambda', newdata=nd)
  expect_equal(dim(pr), c(2,4))

  res <- residuals(fm_C)
  expect_equal(dim(res), dim(y))

  expect_warning(r <- ranef(fm_C))
  expect_equal(dim(r@post), c(4,24,1))

  expect_warning(s <- simulate(fm_C, 2))
  expect_equal(length(s), 2)

  expect_warning(pb <- parboot(fm_C, nsim=1))
  expect_is(pb, "parboot")

  expect_error(gmultmix(~(1|dummy),~1,~1,umf))

})

test_that("gmultmix double model works",{
  # Simulate independent double observer data
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

  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)
  expect_warning(umf <- unmarkedFrameGMM(y=y, obsCovs=list(observer=observer),
         type="double",numPrimary=1))

  # check subset
  umf2 <- umf[1:5,]
  expect_equal(numSites(umf2), 5)

  fm <- gmultmix(~1,~1,~observer, umf)
  expect_equivalent(coef(fm), c(2.2586,0.17385,-0.7425), tol=1e-4)

  gp <- getP(fm)
  expect_equal(dim(gp), c(nSites, 3))

})

test_that("gmultmix dependent double model works",{
  # Simulate independent double observer data
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

  # Fit model
  observer <- matrix(c('A','B'), nSites, 2, byrow=TRUE)

  #expect_error(umf <- unmarkedFrameGMM(y=y, obsCovs=list(observer=observer),
  #       type="depDouble",numPrimary=1))

  expect_warning(umf <- unmarkedFrameGMM(y=y[,1:2], obsCovs=list(observer=observer),
         type="depDouble",numPrimary=1))

  fm <- gmultmix(~1,~1,~observer, umf)
  expect_equivalent(coef(fm), c(1.7762,0.2493,0.2008), tol=1e-4)

  gp <- getP(fm)
  expect_equal(dim(gp), c(nSites, 2))

})

test_that("MRR custom piFun works",{

  alfl <- read.csv(system.file("csv", "alfl.csv", package="unmarked"))
  alfl.covs <- read.csv(system.file("csv", "alflCovs.csv",package="unmarked"),
  row.names=1)
  alfl$captureHistory <- paste(alfl$interval1, alfl$interval2, alfl$interval3,
  sep="")
  alfl$captureHistory <- factor(alfl$captureHistory,
     levels=c("001", "010", "011", "100", "101", "110", "111"))
  alfl$id <- factor(alfl$id, levels=rownames(alfl.covs))

  alfl.v1 <- alfl[alfl$survey==1,]
  alfl.H1 <- table(alfl.v1$id, alfl.v1$captureHistory)
  alfl.v2 <- alfl[alfl$survey==2,]
  alfl.H2 <- table(alfl.v2$id, alfl.v2$captureHistory)
  alfl.v3 <- alfl[alfl$survey==3,]
  alfl.H3 <- table(alfl.v3$id, alfl.v3$captureHistory)


  Y<- array(NA, c(50, 3, 7))
  Y[1:50,1,1:7]<- alfl.H1
  Y[1:50,2,1:7]<- alfl.H2
  Y[1:50,3,1:7]<- alfl.H3

  crPiFun <- function(p) {
    p1 <- p[,1]
    p2 <- p[,2]
    p3 <- p[,3]
    cbind("001" = (1 - p1) * (1 - p2) *      p3,
          "010" = (1 - p1) *      p2  * (1 - p3),
          "011" = (1 - p1) *      p2  *      p3,
          "100" =      p1  * (1 - p2) * (1 - p3),
          "101" =      p1  * (1 - p2) *      p3,
          "110" =      p1  *      p2  * (1 - p3),
          "111" =      p1  *      p2  *      p3)
  }

  intervalMat <- matrix(c('1','2','3'), 50, 3, byrow=TRUE)
  class(alfl.H1) <- "matrix"
  o2y <- matrix(1, 3, 7)

  ywide<- as.matrix( cbind(alfl.H1, alfl.H2)  )
  umf.cr1 <- unmarkedFrameGMM(y=ywide,
    obsCovs=NULL, yearlySiteCovs=NULL,
      obsToY=o2y, numPrimary=2, piFun="crPiFun")

  expect_equal(dim(umf.cr1@obsToY)[1] , 6)
  expect_equal(dim(umf.cr1@obsToY)[2] , 14)
})

test_that("R and C++ engines give identical results",{
  y <- matrix(0:3, 5, 4)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
  yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))

  umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
  expect_warning(fm_R <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23,
                                  engine="R", control=list(maxit=1)))
  expect_warning(fm_C <- gmultmix(~x, ~yr, ~o1 + o2, data = umf, K=23,
                                  engine="C", control=list(maxit=1)))
  expect_equal(coef(fm_R), coef(fm_C))


})
