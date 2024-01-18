context("distsampOpen fitting function")
skip_on_cran()

simData <- function(lambda=1, gamma=0.5, omega=0.8, sigma=40, scale=NULL,
                    M=100, T=5, J=4,type="line", keyfun="halfnorm")
{
    y <- array(NA, c(M, J, T))
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
if(keyfun=="halfnorm"){
    if(type=="point")
       g <- function(x, sig) exp(-x^2/(2*sig^2))*x
    if(type=="line")
       g <- function(x, sig) exp(-x^2/(2*sig^2))
}
if(keyfun=="exp"){
    if(type=="point")
       g <- function(x, sig) exp(-x/sig)*x
    if(type=="line")
       g <- function(x, sig) exp(-x/sig)
}
if(keyfun=="hazard"){
    shape<- sigma
    if(type=="point")
       g <- function(x, shape, scale) (1 - exp(-(r/shape)^-scale)) * x
    if(type=="line")
       g <- function(x, sig) (1 - exp(-(x/shape)^-scale) )
}
if(keyfun=="uniform"){
     if(type=="point")
       g <-  function(x,sig) x
    if(type=="line")
       g <- function(x,sig) 1
}

    cp <- u <- a <- numeric(J)
if(keyfun!="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]
    cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value
    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value
    }
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

}
if(keyfun=="uniform"){
if(type=="point"){
   a[1] <- pi*db[2]^2
   for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])

    }
    }
if(type=="line"){
   L <-  1
   a[1] <- L*db[2]

    for(j in 2:J) {
      a[j] <-  db[j+1]  - sum(a[1:j])

    }
    }
    u <- a / sum(a)
    cp <- a
    cp[j+1] <- 0

}
    for(i in 1:M) {
    N[i,1] <- rpois(1, lambda)

    y[i,1:J,1] <- rmultinom(1, N[i,1], cp)[1:J]

    for(t in 1:(T-1)) {
        S[i,t] <- rbinom(1, N[i,t], omega)
        G[i,t] <- rpois(1, gamma)
        N[i,t+1] <- S[i,t] + G[i,t]
        y[i,1:J,t+1] <- rmultinom(1, N[i,t+1], cp)[1:J]
        }
    }
    #cp <- array(cp, c(J, M, T))
    #cp <- matrix(aperm(cp, c(2,1,3)), M)
    return(list(y=matrix(y, M),N=N,cp=cp))
}

test_that("unmarkedFrameDSO build properly", {
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=15,type="line",
            keyfun="halfnorm")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=15,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))
  expect_error(unmarkedFrameDSO(y=y,numPrimary=15,
                  dist.breaks=c(0,25,50,75,100), survey="line", tlength=c(1,1)))

  expect_error(unmarkedFrameDSO(y=y,numPrimary=15,
                  dist.breaks=c(0,25,50,75,100), survey="point",
                  tlength=rep(1,100)))
  expect_error(unmarkedFrameDSO(y=y, numPrimary=15, dist.breaks=c(25,50,75,100),
                                  survey='line', unitsIn='m', tlength=rep(1,100)))

  # subset sites
  umf_sub <- umf[1:3,]
  expect_equal(nrow(umf_sub@y), 3)
})

test_that("dso halfnorm key function works",{
  set.seed(456)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=50, T=7,type="line",
            keyfun="halfnorm")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=7,
                          siteCovs=data.frame(x1=rnorm(50)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 50))

  # Check K that is too small
  expect_error(distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=5,keyfun="halfnorm"))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=10,keyfun="halfnorm")

  expect_equivalent(coef(fm), c(1.4185,1.0471,-0.8275,3.1969,-0.0790),
                     tol=1e-4)

  pr <- predict(fm, type='lambda')
  expect_equal(as.numeric(pr[1,]),
                     c(4.1308,0.4965,3.2622,5.2306), tol=1e-4)

  pval <- getP(fm)
  expect_equal(dim(pval), dim(y))
  expect_equal(pval[1,1:4], c(0.2110,0.0776,0.0104,0.0005),
                     tol=1e-4)

  r <- residuals(fm)
  expect_equal(dim(r), dim(y))
  expect_equal(r[1,1:2], c(-0.8717,-0.3207),tol=1e-4)

  ran <- ranef(fm)
  expect_equal(bup(ran)[1,1], 2.8916, tol=1e-4)

  set.seed(123)
  sim <- simulate(fm, nsim=2)
  expect_equal(length(sim), 2)
  expect_equal(sum(sim[[1]]), 442)

  fm2 <- update(fm, pformula=~1)
  expect_equal(length(coef(fm2)), 4)

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=10,keyfun="halfnorm",
                     mixture="ZIP")
  expect_warning(pr <- predict(fm, 'lambda'))
  expect_is(pr, "data.frame")
  expect_equal(as.numeric(pr[1,]), c(4.0742, NA,NA,NA), tol=1e-4)

  #Point
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=50, T=7,type="point",
            keyfun="halfnorm")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=7,
                          siteCovs=data.frame(x1=rnorm(50)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="point",
                          unitsIn="m")

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=10,keyfun="halfnorm")
  expect_equivalent(coef(fm), c(1.2766,0.4651,0.23799,3.24919,0.00821),
                                 tol=1e-4)

  #Check error with random effects in formula
  expect_error(distsampOpen(~(1|dummy), ~1, ~1, ~1, data=umf, K=30))
})

test_that("distsampOpen works with NAs", {
  set.seed(456)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=15, T=2,type="line",
            keyfun="halfnorm")$y

  y[5,1:4] <- NA
  y[2,1] <- NA
  sc = data.frame(x1=rnorm(15))
  sc$x1[3] <- NA
  ysc = data.frame(x2=rnorm(15*2))
  ysc$x2[1] <- NA

  umf <- unmarkedFrameDSO(y = y, numPrimary=2,
            siteCovs=sc,
            yearlySiteCovs=ysc,
            dist.breaks = c(0, 25, 50, 75, 100), survey="line",
            unitsIn="m",tlength=rep(1, 15))

  expect_warning(fm <- distsampOpen(~x1, ~x2, ~1, ~1, data=umf, K=7, keyfun="halfnorm"))
  expect_equivalent(coef(fm), c(1.3058,-0.2966,-7.9133,-7.9281,8.6582,3.3108), tol=1e-4)

  set.seed(123)
  ysim <- simData(lambda=5, gamma=2, omega=0.5, sigma=40, M=50, T=5,type="line",
            keyfun="halfnorm")
  y <- ysim$y
  y[2,1] <- NA

  umf <- unmarkedFrameDSO(y = y, numPrimary=5,
            dist.breaks = c(0, 25, 50, 75, 100), survey="line",
            unitsIn="m",tlength=rep(1, 50))

  fm <- distsampOpen(~1, ~1, ~1, ~1, data=umf, K=10, keyfun="halfnorm")

  expect_warning(r <- ranef(fm))
  expect_equal(cor(bup(r)[,1],ysim$N[,1]), 0.6593, tol=1e-4)

})

test_that("distsampOpen exp keyfunction works", {
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=50, T=5,type="line",
            keyfun="exp")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=5,
                          siteCovs=data.frame(x1=rnorm(50)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 50))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=10,keyfun="exp")
  expect_equivalent(coef(fm), c(1.47976,0.38259,0.15922,3.2837,-0.01028),
                     tol=1e-4)

  #Point
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=30, T=3,type="point",
            keyfun="exp")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=3,
                          siteCovs=data.frame(x1=rnorm(30)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="point",
                          unitsIn="m")

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=10,keyfun="exp")
  expect_equivalent(coef(fm), c(1.5032,-3.6228,9.25622,3.09501,0.05183),
                     tol=1e-4)

})

test_that("distsampOpen uniform keyfun works", {
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=50, T=5,type="line",
            keyfun="uniform")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=5,
                          siteCovs=data.frame(x1=rnorm(50)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 50))

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=15,keyfun="unif")
  expect_equivalent(coef(fm), c(1.4586,0.7262,-0.05239),
                     tol=1e-4)

  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=50, T=5,type="point",
            keyfun="uniform")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=5,
                          siteCovs=data.frame(x1=rnorm(50)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="point",
                          unitsIn="m")

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=15,keyfun="unif")
  expect_equivalent(coef(fm), c(1.4730,0.6887,-0.1756),
                     tol=1e-4)
})

test_that("distsampOpen hazard keyfun works", {
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, scale=1,
               M=30, T=3,type="line", keyfun="hazard")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=3,
                          siteCovs=data.frame(x1=rnorm(30)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 30))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=30,keyfun="hazard")
  expect_equivalent(coef(fm), c(1.2671,0.6878,-1.715,4.0409,0.00578,-0.08012),
                     tol=1e-4)
})

test_that("distsampOpen NB abundance model works",{

  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=30, T=3,type="line",
            keyfun="exp")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=3,
                          siteCovs=data.frame(x1=rnorm(30)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 30))

  fm <- distsampOpen(~1, ~1, ~1, ~x1, data = umf, K=10,keyfun="exp",
                 mixture="NB")
  expect_equivalent(coef(fm), c(1.0517,0.5117,0.01159,3.34406,0.14606,8.23921), tol=1e-4)

})

test_that("distsampOpen dynamics models work",{
  set.seed(123)
  y <- simData(lambda=4, gamma=2, omega=0.5, sigma=25, M=100, T=10,type="line",
            keyfun="uniform")$y
  umf <- unmarkedFrameDSO(y = y, numPrimary=10,
                          siteCovs=data.frame(x1=rnorm(100)),
                          dist.breaks = c(0, 25, 50, 75, 100), survey="line",
                          unitsIn="m",tlength=rep(1, 100))

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25,keyfun="unif",
                    dynamics="notrend")
  expect_equivalent(coef(fm), c(1.4080889, -0.1006024), tol=1e-5)

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="trend")
  expect_equivalent(coef(fm), c(1.518695, -0.0143889), tol=1e-5)

  fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="autoreg")
  expect_equivalent(coef(fm), c(1.518686, -0.018026, -5.628779), tol=1e-4)

  #Sketchy estimates
  #Maybe just because data were simulated using a different process?
  #Leaving these in for now just to make sure they run without errors
  expect_warning(fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="gompertz"))

  expect_warning(fm <- distsampOpen(~1, ~1, ~1, data = umf, K=25, keyfun="unif",
                     dynamics="ricker"))

})
