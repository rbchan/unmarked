context("gdistremoval fitting function")
skip_on_cran()

simData <- function(lambda=1, sigma=40, scale=NULL, remP=0.4, remJ=4,
                    M=100, J=4,type="point", keyfun="halfnorm", T=1, phi=1)
{
    y <- array(NA, c(M, T*J))
    yRem <- matrix(NA, M, T*remJ)
    N <- rep(NA, M)

piRem <- rep(NA,remJ+1)
piRem[1] <- remP
for (j in 2:remJ){
  piRem[j] <- piRem[j-1] / remP * (1-remP) * remP
}
piRem[remJ+1] <- 1-sum(piRem, na.rm=T)

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

      N[i] <- rpois(1, lambda)

      yRem_sub <- y_sub <- c()

      for (t in 1:T){

        yRemT <- rep(0, remJ)
        yT <- rep(0, J)

        dist_class <- rem_class <- rep(NA, N[i])

        prob_in <- sum(piRem[1:remJ])*sum(cp[1:J])*phi

        if(N[i] > 0){
        for (n in 1:N[i]){

          keep <- rbinom(1,1,prob_in)
          if(keep==1){
            rem_class[n] <- sample(1:remJ, 1, prob=piRem[1:remJ]/sum(piRem[1:remJ]))
            dist_class[n] <- sample(1:J, 1, prob=cp[1:J]/sum(cp[1:J]))

            yT[dist_class[n]] <- yT[dist_class[n]] + 1
            yRemT[rem_class[n]] <- yRemT[rem_class[n]] + 1
          }
        }
        }
        yRem_sub <- c(yRem_sub, yRemT)
        y_sub <- c(y_sub, yT)
      }
      yRem[i,] <- yRem_sub
      y[i,] <- y_sub
    }

    return(list(y=y,N=N,yRem=yRem))
}

simDataRand <- function(lambda=1, lam_sd=1, groups=10,
                        sigma=40, remP=0.4, remJ=4,
                        M=100, J=4, T=1, phi=1)
{
    y <- array(NA, c(M, T*J))
    yRem <- matrix(NA, M, T*remJ)
    N <- rep(NA, M)

group_eff <- rnorm(groups, 0, lam_sd)
names(group_eff) <- letters[1:groups]
group <- sample(letters[1:groups], M, replace=T)
rand_ef <- group_eff[group]
lambda <- exp(log(lambda) + rand_ef)

piRem <- rep(NA,remJ+1)
piRem[1] <- remP
for (j in 2:remJ){
  piRem[j] <- piRem[j-1] / remP * (1-remP) * remP
}
piRem[remJ+1] <- 1-sum(piRem, na.rm=T)

    db <- c(0, 25, 50, 75, 100)
    if(length(db)-1 != J)
      stop("hey, what")
g <- function(x, sig) exp(-x^2/(2*sig^2))*x

cp <- u <- a <- numeric(J)

  a[1] <- pi*db[2]^2
  cp[1] <- integrate(g, db[1], db[2], sig=sigma)$value * 2 * pi
    for(j in 2:J) {
      a[j] <- pi*db[j+1]^2 - sum(a[1:j])
      cp[j] <- integrate(g, db[j], db[j+1], sig=sigma)$value * 2*pi
    }
    u <- a / sum(a)
    cp <- cp / a * u
    cp[j+1] <- 1-sum(cp)

    for(i in 1:M) {

      N[i] <- rpois(1, lambda[i])

      yRem_sub <- y_sub <- c()

      for (t in 1:T){

        yRemT <- rep(0, remJ)
        yT <- rep(0, J)

        dist_class <- rem_class <- rep(NA, N[i])

        prob_in <- sum(piRem[1:remJ])*sum(cp[1:J])*phi

        if(N[i] > 0){
        for (n in 1:N[i]){

          keep <- rbinom(1,1,prob_in)
          if(keep==1){
            rem_class[n] <- sample(1:remJ, 1, prob=piRem[1:remJ]/sum(piRem[1:remJ]))
            dist_class[n] <- sample(1:J, 1, prob=cp[1:J]/sum(cp[1:J]))

            yT[dist_class[n]] <- yT[dist_class[n]] + 1
            yRemT[rem_class[n]] <- yRemT[rem_class[n]] + 1
          }
        }
        }
        yRem_sub <- c(yRem_sub, yRemT)
        y_sub <- c(y_sub, yT)
      }
      yRem[i,] <- yRem_sub
      y[i,] <- y_sub
    }

    return(list(y=y,N=N,yRem=yRem,group=group,group_eff=group_eff))
}

test_that("unmarkedFrameGDR is constructed correctly",{
  set.seed(123)

  # Single primary period
  simdat <- simData(remJ=5)
  sc <- data.frame(sc1=rnorm(nrow(simdat$y)))
  oc <- data.frame(oc1=rnorm(length(simdat$yRem)))

  umf <- unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                          unitsIn='m')
  expect_true(inherits(umf, "unmarkedFrameGDR"))

  # Check subsetting
  umf_sub <- umf[1:3,]
  expect_true(inherits(umf_sub, "unmarkedFrameGDR"))
  expect_equivalent(numSites(umf_sub), 3)
  expect_error(umf[,1:2])

  # Input mistake handling

  # Wrong number of dist.breaks
  expect_error(unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30),
                          unitsIn='m'))

  # Wrong number of period.lengths
  expect_error(unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                          unitsIn='m', period.lengths=c(1,1,1)))

  # row sums of yDistance and yRemoval don't match
  yRem_bad <- simdat$yRem
  yRem_bad[which(yRem_bad>0)[1]] <- 0
  expect_error(unmarkedFrameGDR(simdat$y, yRem_bad, siteCovs=sc,
                          obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                          unitsIn='m'))

  # Multi primary period
  simdat <- simData(T=3, remJ=5)
  sc <- data.frame(sc1=rnorm(nrow(simdat$y)))
  oc <- data.frame(oc1=rnorm(length(simdat$yRem)))
  ysc <- data.frame(ysc1=rnorm(100*3))

  umf2 <- unmarkedFrameGDR(simdat$y, simdat$yRem, siteCovs=sc,
                           obsCovs=oc, dist.breaks=c(0,10,20,30,40),
                           yearlySiteCovs=ysc, unitsIn='m', numPrimary=3)

  # Check subsetting
  umf_sub <- umf2[1:3,]
  expect_true(inherits(umf_sub, "unmarkedFrameGDR"))
  expect_equivalent(numSites(umf_sub), 3)
  expect_equivalent(nrow(siteCovs(umf_sub)), 3)
  expect_equivalent(nrow(obsCovs(umf_sub)), 3*15)
  expect_equivalent(nrow(yearlySiteCovs(umf_sub)), 3*3)
  umf_sub <- umf2[,1:2]
  expect_equivalent(ncol(umf_sub@yRemoval), 10)
  expect_equivalent(ncol(umf_sub@yDistance), 8)
  expect_equivalent(nrow(umf_sub@obsCovs), 10*100)
  expect_equivalent(nrow(umf_sub@yearlySiteCovs), 2*100)
  expect_equivalent(umf_sub@numPrimary, 2)
})

test_that("gdistremoval can fit models",{
  set.seed(123)

  sc <- data.frame(sc1=rnorm(50))
  oc <- data.frame(oc1=rnorm(5*50))

  # Half-normal
  dat <- simData(lambda=5, sigma=50, M=50, J=4, remP=0.2, remJ=5)
  umf <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')

  fit <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf)
  expect_is(fit, "unmarkedFitGDR")
  expect_equivalent(coef(fit), c(1.4571,0.3374,4.0404,-1.65389,0.16789), tol=1e-3)

  # With unequal period lengths
  umf2 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m',
                         period.lengths=c(1,5,1,1,1))

  fit2 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf2)
  expect_equivalent(coef(fit2), c(2.7732,0.3241,4.0477,-3.9774,0.1418), tol=1e-3)

  # With negative binomial
  fit3 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf,
                       mixture="NB")
  expect_equivalent(coef(fit3), c(1.4571,0.33742,8.3535,4.0404,-1.6539,0.1679), tol=1e-3)

  # With exponential
  set.seed(123)
  dat <- simData(lambda=5, sigma=50, M=50, J=4, remP=0.2, remJ=5, keyfun="exp")
  umf4 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit4 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf4,
                       keyfun="exp")
  expect_equivalent(coef(fit4), c(1.54527,0.0045,4.2135,-1.8776,-0.27996), tol=1e-3)

  # With hazard
  set.seed(123)
  dat <- simData(lambda=5, sigma=50, M=50, J=4, remP=0.2, remJ=5)
  umf5 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit5 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf5,
                       keyfun="hazard")
  expect_equivalent(coef(fit5), c(1.4477,0.0569,4.10349,1.39400,-1.419985,-0.08351), tol=1e-3)

  # With uniform
  set.seed(123)
  dat <- simData(lambda=5, sigma=50, M=50, J=4, remP=0.2, remJ=5)
  umf6 <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m')
  fit6 <- gdistremoval(~sc1,removalformula=~oc1,distanceformula=~1, data=umf6,
                       keyfun="uniform")
  expect_equivalent(coef(fit6), c(0.7887,0.0569,-1.4197,-0.083708), tol=1e-3)

  # Methods
  gp <- getP(fit)
  expect_equivalent(dim(gp$dist), c(50,4,1))
  expect_equivalent(dim(gp$rem), c(50,5,1))

  s <- simulate(fit, 2)
  expect_equivalent(length(s), 2)
  expect_equivalent(dim(s[[1]]$yDistance), dim(fit@data@yDistance))
  expect_equivalent(dim(s[[1]]$yRemoval), dim(fit@data@yRemoval))

  r <- ranef(fit)
  expect_equivalent(length(bup(r)), 50)

  pb <- parboot(fit, nsim=2)
  expect_is(pb, "parboot")

  # Fit list construction
  fl <- fitList(fits=list(fit1=fit, fit2=fit))
  expect_is(fl, "unmarkedFitList")
  ms <- modSel(fl)
  expect_is(ms, "unmarkedModSel")

  #plotEffects
  pdf(NULL)
  plotEffects(fit, 'lambda', "sc1")
  plotEffects(fit, 'rem', 'oc1')
  dev.off()
})

test_that("gdistremoval predict method works",{

  set.seed(123)

  sc <- data.frame(sc1=rnorm(50), sc2=sample(letters[1:5],50,replace=T))
  oc <- data.frame(oc1=rnorm(5*50))

  # Half-normal
  dat <- simData(lambda=5, sigma=50, M=50, J=4, remP=0.2, remJ=5)
  umf <- expect_warning(unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m'))

  fit <- gdistremoval(~sc1+sc2,removalformula=~oc1+sc1,distanceformula=~1, data=umf)

  pr <- predict(fit, 'lambda')
  expect_equivalent(dim(pr), c(50,4))
  nd <- data.frame(sc1=c(0,1), sc2='a')
  pr <- predict(fit, 'lambda', newdata=nd)
  expect_equivalent(dim(pr), c(2,4))

  nd <- data.frame(sc1=c(0,1), sc2=letters[6])
  expect_error(predict(fit, 'lambda', newdata=nd))

  pr <- predict(fit, 'rem')
  expect_equivalent(dim(pr), c(5*50,4))
  nd <- data.frame(oc1=c(0,1))
  expect_error(predict(fit, 'lambda', newdata=nd))
  nd <- data.frame(oc1=c(0,1), sc1=c(0,1))
  pr <- predict(fit, 'rem', newdata=nd)
  expect_equivalent(dim(pr), c(2,4))

  pr <- predict(fit, 'dist')
  expect_equivalent(dim(pr), c(50,4))

  expect_error(predict(fit, 'test'))
})

test_that("gdistremoval handles NAs",{
  set.seed(123)

  sc <- data.frame(sc1=rnorm(50), sc2=sample(letters[1:5],50,replace=T))
  oc <- data.frame(oc1=rnorm(5*50))

  # Half-normal
  dat <- simData(lambda=5, sigma=50, M=50, J=4, remP=0.2, remJ=5)

  yDist <- dat$y
  yDist[1,1] <- NA
  yDist[2,] <- NA

  yRem <- dat$yRem
  yRem[3,1] <- NA
  yRem[2,] <- NA

  expect_warning(umf <- unmarkedFrameGDR(yDist, yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m'))

  fit <- gdistremoval(~1,removalformula=~1,distanceformula=~1, data=umf)
  expect_equivalent(coef(fit), c(2.0675,3.908,-2.1433), tol=1e-3)

  # Can't have missing site covs
  umf2 <- umf
  umf2@siteCovs$sc1[1] <- NA
  expect_error(gdistremoval(~sc1,removalformula=~1,distanceformula=~1, data=umf2))
  
  # This errors because missing obs cov does not match missing removal data
  umf2 <- umf
  umf2@obsCovs$oc1[1] <- NA
  expect_error(gdistremoval(~1,removalformula=~oc1,distanceformula=~1, data=umf2))

  # This does not error because missing obs cov matches missing removal data
  umf2 <- umf
  umf2@obsCovs$oc1[6] <- NA
  fit <- gdistremoval(~1,removalformula=~oc1,distanceformula=~1, data=umf2)
  expect_true(is.na(predict(fit, 'rem')$Predicted[6]))
})

test_that("multi-period data works with gdistremoval",{
  set.seed(123)

  sc <- data.frame(sc1=rnorm(30))
  oc <- data.frame(oc1=rnorm(5*30*5))
  ysc <- data.frame(ysc1=rnorm(30*5))

  dat <- simData(lambda=5, sigma=30, M=30, T=5, J=4, remP=0.2, remJ=5, phi=0.5)
  umf <- unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc, yearlySiteCovs=ysc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m', numPrimary=5)

  fit <- gdistremoval(~sc1,phiformula=~ysc1, removalformula=~oc1,distanceformula=~1, data=umf)

  expect_equivalent(coef(fit),
                     c(2.1013,-0.1142,-1.3187,-0.1483,3.3981,-0.5142,0.233678),
                     tol=1e-3)

  # Predict
  pr <- predict(fit, 'phi')
  expect_equal(dim(pr), c(30*5,4))
  expect_equal(as.numeric(pr[1,1]), c(0.1916), tol=1e-3)

  # getP
  gp <- getP(fit)
  expect_equal(lapply(gp, dim)$phi, c(30,5))

  # ranef
  r <- ranef(fit)
  expect_equal(dim(r@post), c(30, 44, 1))
  expect_equal(length(bup(r)), 30)

  # fitted
  ft <- fitted(fit)
  expect_equal(dim(ft$dist), dim(fit@data@yDistance))

  # Entire missing secondary period
  umf2 <- umf
  # remove 2nd period at first site
  umf2@yDistance[1,5:8] <- NA
  umf2@yRemoval[1, 6:10] <- NA
  umf2@obsCovs$oc1[6:10] <- NA
  umf2@yearlySiteCovs$ysc1[2] <- NA
  

  fit2 <- gdistremoval(~sc1,phiformula=~ysc1, removalformula=~oc1,distanceformula=~1, data=umf2)

  gp <- getP(fit2)
  expect_true(is.na(gp$phi[1,2]))

  pr <- predict(fit2, 'rem', level=NULL)
  expect_true(all(is.na(pr$Predicted[6:10])))

  r2 <- ranef(fit2)
  expect_true(!any(is.na(r2@post)))

  s <- simulate(fit2)
  expect_true(all(is.na(s[[1]]$yDistance[1,5:8])))
  
  res <- residuals(fit2)
  expect_true(all(is.na(res$rem[1,6:10])))
  pb <- parboot(fit2, nsim=2)
  expect_is(pb, "parboot")
})

test_that("gdistremoval works with random effects",{

  set.seed(123)
  dat <- simDataRand(lambda=5, lam_sd=0.4, groups=10, sigma=50, M=50, J=4, remP=0.2, remJ=5) #
  sc <- data.frame(sc1=rnorm(50), group=dat$group)
  oc <- data.frame(oc1=rnorm(5*50))
  umf <- expect_warning(unmarkedFrameGDR(dat$y, dat$yRem, siteCovs=sc, obsCovs=oc,
                         dist.breaks=c(0,25,50,75,100), unitsIn='m'))
  fit <- gdistremoval(~sc1+(1|group),removalformula=~oc1,distanceformula=~1, data=umf)

  expect_equivalent(coef(fit), c(1.50680,0.07968,3.9373,-1.2148,0.09809),
                     tol=1e-4)
  expect_equal(sigma(fit)$sigma, 0.2866, tol=1e-4)

  pr <- predict(fit, "lambda")
  expect_is(pr, "data.frame")
  pr <- predict(fit, "lambda", newdata=umf@siteCovs[1:2,])
  expect_is(pr, "data.frame")

  s <- simulate(fit, 2)
  expect_is(s, "list")

  pb <- parboot(fit, nsim=1)
  expect_is(pb, "parboot")

})
