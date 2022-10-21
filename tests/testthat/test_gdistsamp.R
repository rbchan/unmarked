context("gdistsamp fitting function")
skip_on_cran()

test_that("unmarkedFrameGDS subset works",{
    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))

    umf1 <- unmarkedFrameGDS(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        numPrimary = 3,
        survey="point",
        dist.breaks=c(0, 10, 20, 30),
        unitsIn="m")

    dat <- as(umf1, "data.frame")
    expect_equal(nrow(dat), nrow(y))

    umf1.site1 <- umf1[1,]
    expect_equal(umf1.site1@y, y[1,, drop=FALSE])
    expect_equal(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    expect_equivalent(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    expect_equal(umf1.site1@numPrimary, 3)
    expect_equal(umf1.site1@survey, "point")

    umf1.sites1and3 <- umf1[c(1,3),]


    umf2 <- unmarkedFrameGDS(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        numPrimary = 3,
        survey="line",
        dist.breaks=c(0, 10, 20, 30),
        tlength=rep(1,nrow(y)),
        unitsIn="m")

    dat <- as(umf2, "data.frame")

    umf2.site1 <- umf2[1,]
    expect_equal(umf2.site1@y, y[1,, drop=FALSE])
    expect_equal(umf2.site1@siteCovs, sc[1,, drop=FALSE])
    expect_equivalent(unlist(umf2.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    expect_equal(umf2.site1@numPrimary, 3)
    expect_equal(umf2.site1@survey, "line")

    umf2.sites1and3 <- umf2[c(1,3),]
})

test_that("gdistsamp with halfnorm keyfunction works",{
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    # R and C give same result
    fm_R <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE, engine="R",
                      control=list(maxit=1))
    fm_C <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE, engine="C",
                      control=list(maxit=1))
    expect_equal(coef(fm_R), coef(fm_C))

    #When output = density
    #fm_R <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE, engine="R")
    fm_C <- gdistsamp(~1, ~1, ~1, umf, output="density", se=FALSE, engine="C")

    expect_equivalent(coef(fm_C),c(0.6584008,-0.4440278,3.4817094),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #When output = abundance
    #fm_R <- gdistsamp(~1, ~1, ~1, umf, output="abund", se=FALSE, engine="R")
    fm_C <- gdistsamp(~1, ~1, ~1, umf, output="abund", se=FALSE, engine="C")

    expect_equivalent(coef(fm_C),c(1.35067,-0.44262,3.48149),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #With covariates
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    #fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE, engine="C")

    expect_equivalent(coef(fm_C),c(1.24510,0.54419,-1.28146,-0.109737,
                                    3.46295,-0.13228),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    # methods
    res <- residuals(fm_C)
    expect_equal(dim(res), c(30,15))
    expect_equal(res[1,1], -0.07133, tol=1e-4)
    r <- ranef(fm_C)
    expect_equal(dim(r@post), c(30,108,1))
    expect_equal(bup(r)[1], 1.94300, tol=1e-5)
    s <- simulate(fm_C, 2)
    expect_is(s, "list")
    expect_equal(length(s), 2)
    pb <- parboot(fm_C, nsim=1)
    expect_is(pb, "parboot")

    #Negative binomial
    #fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
    #                  mixture="NB", engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      mixture="NB", engine="C")

    expect_equivalent(coef(fm_C),c(1.41241,0.52442,-1.49024,-0.10546,
                                    3.46284,-0.129831,3.16892),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #With missing values
    yna <- y
    yna[1,c(1,6)] <- NA
    umf <- unmarkedFrameGDS(y = yna, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    #fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
    #                  se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
                      se=FALSE, engine="C")

    expect_equivalent(coef(fm_C),c(1.35065,0.52558,-1.39758,-0.10675,
                                    3.46283,-0.136344),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #With an entire session missing
    yna <- y
    yna[1,1:J] <- NA
    umf <- unmarkedFrameGDS(y = yna, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    #fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
    #                  se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density",
                      se=FALSE, engine="C")

    expect_equivalent(coef(fm_C),c(1.30815,0.53527,-1.35387,-0.11038,
                                    3.46293,-0.13458),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point
    set.seed(123)
    data(issj)
    covs <- issj[1:100,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[1:100,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    #fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
    #                  engine="R", se=F)
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      engine="C", se=F)

    expect_equivalent(coef(fm_C),c(-4.1009,-0.2363,3.5557,8.7543),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    # Check error when formula has random effects
    expect_error(gdistsamp(~(1|dummy), ~1, ~1, umf))

    # R and C engines return same result
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      engine="C", se=F, control=list(maxit=1))
    fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      engine="R", se=F, control=list(maxit=1))
    expect_equal(coef(fm_C), coef(fm_R))


})

test_that("gdistsamp with uniform keyfunction works",{
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    # R and C engines return same result
    fm_R <- gdistsamp(~par1, ~par2, ~1, umf, output="density",
                      keyfun="uniform", se=FALSE, engine="C", control=list(maxit=1))
    fm_C <- gdistsamp(~par1, ~par2, ~1, umf, output="density",
                      keyfun="uniform", se=FALSE, engine="R", control=list(maxit=1))
    expect_equal(coef(fm_R), coef(fm_C))

    #fm_R <- gdistsamp(~par1, ~par2, ~1, umf, output="density",
    #                  keyfun="uniform", se=FALSE, engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~1, umf, output="density",
                      keyfun="uniform", se=FALSE, engine="C")

    expect_equivalent(coef(fm_C),c(1.17120,0.54748,-1.60963,-0.13009),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point: doesn't work with this dataset, find another one?
    #OR maybe uniform just doesn't work with point?
    set.seed(123)
    data(issj)
    covs <- issj[,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    expect_error(gdistsamp(~1, ~1, ~1, jayumf, output='density',
                      keyfun='uniform', engine="R", se=F))
    expect_error(gdistsamp(~elevation, ~1, ~1, jayumf, output='density',
                      keyfun='uniform', engine="C", se=F))
})


test_that("gdistsamp with exp keyfunction works",{
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    #With covariates
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    # R and C engines return same result
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="exp",engine="C", control=list(maxit=1))
    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="exp",engine="R", control=list(maxit=1))
    expect_equal(fm_C@AIC, fm_R@AIC)

    #fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
    #                  keyfun="exp",engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="exp",engine="C")

    expect_equivalent(coef(fm_C),c(1.28243,0.54312,-1.16608,-0.101122,
                                    3.86666,-0.2492846),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point
    set.seed(123)
    data(issj)
    covs <- issj[1:100,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[1:100,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    #fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
    #                  keyfun="exp",engine="R", se=F)
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="exp",engine="C", se=F)

    expect_equivalent(coef(fm_C),c(-3.1952,-0.1244,3.7986,3.6574),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="exp",engine="C", se=F,control=list(maxit=1))
    fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="exp",engine="R", se=F, control=list(maxit=1))
    expect_equal(fm_C@AIC, fm_R@AIC, tol=1e-5)

})

test_that("gdistsamp with hazard keyfunction works",{
    #Line
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    #With covariates
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    # R and C engines give same result
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="hazard",engine="C", control=list(maxit=1))
    fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="hazard",engine="R", control=list(maxit=1))
    expect_equal(fm_C@AIC, fm_R@AIC, tol=1e-5)

    #fm_R <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
    #                  keyfun="hazard",engine="R")
    fm_C <- gdistsamp(~par1, ~par2, ~par3, umf, output="density", se=FALSE,
                      keyfun="hazard",engine="C")

    expect_equivalent(coef(fm_C),c(1.29425,0.54233,-1.41658,-0.09267,
                                    3.45436,-0.19978,0.8270215),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-4)

    #Point
    set.seed(123)
    data(issj)
    covs <- issj[1:100,c("elevation","forest","chaparral")]
    area <- pi*300^2 / 100^2
    # Area in ha
    jayumf <- unmarkedFrameGDS(y=as.matrix(issj[1:100,1:3]),
                               siteCovs=data.frame(covs, area),
                               yearlySiteCovs=data.frame(covs),
                               numPrimary=1,
                              dist.breaks=c(0, 100, 200, 300), unitsIn="m",
                              survey="point")
    sc <- siteCovs(jayumf)
    sc.s <- scale(sc)
    sc.s[,"area"] <- pi*300^2 / 10000 # Don't standardize area
    covs <- siteCovs(jayumf) <- sc.s
    #fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
    #                  keyfun="hazard",engine="R", se=F)
    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="hazard",engine="C", se=F)

    expect_equivalent(coef(fm_C),c(0.8584,-0.04336,-0.70738,3.2762,0.1807),tol=1e-4)
    #expect_equivalent(coef(fm_R),coef(fm_C),tol=1e-3)

    fm_C <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="hazard",engine="C", se=F,control=list(maxit=1))
    fm_R <- gdistsamp(~elevation, ~1, ~chaparral, jayumf, output='density',
                      keyfun="hazard",engine="R", se=F, control=list(maxit=1))
    expect_equal(fm_C@AIC, fm_R@AIC, tol=1e-3)
})

test_that("predict works for gdistsamp",{
    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    ysc <- data.frame(fac_cov = factor(rep(letters[1:T], R)))

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=ysc,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    fm_C <- gdistsamp(~1, ~fac_cov -1, ~1, umf,
                      keyfun="halfnorm", se=TRUE, engine="C")

    #lambda
    pr <- predict(fm_C, "lambda")

    expect_equivalent(dim(pr), c(30,4))
    expect_equivalent(pr[1,1], 3.767935, tol=1e-5)
    nd <- data.frame(par1=0)
    pr2 <- predict(fm_C, type='lambda', newdata=nd)

    expect_equivalent(dim(pr2), c(1,4))
    expect_equivalent(pr2[1,1], 3.767935, tol=1e-5)

    #phi
    pr <- predict(fm_C, "phi")

    expect_equivalent(dim(pr), c(90,4))
    expect_equivalent(pr[1,1], 0.4461197, tol=1e-5)

    nd <- data.frame(fac_cov=factor(letters[1:3]))
    pr2 <- predict(fm_C, type="phi", newdata=nd)

    expect_equivalent(dim(pr2), c(3,4))
    expect_equivalent(pr2[1,1], 0.4461197, tol=1e-5)

    #sigma
    pr <- predict(fm_C, "det")

    expect_equivalent(dim(pr), c(90,4))
    expect_equivalent(pr[1,1], 32.51537, tol=1e-5)

    nd <- data.frame(par3=0)
    pr2 <- predict(fm_C, type='det', newdata=nd)

    expect_equivalent(dim(pr2), c(1,4))
    expect_equivalent(pr2[1,1], 32.51537, tol=1e-5)
})

test_that("gdistsamp handles NAs",{
  set.seed(343)
  R <- 30 # number of transects
  T <- 3  # number of replicates
  strip.width <- 50
  transect.length <- 100
  breaks <- seq(0, 50, by=10)

  lambda <- 5 # Abundance
  phi <- 0.6  # Availability
  sigma <- 30 # Half-normal shape parameter

  J <- length(breaks)-1
  y <- array(0, c(R, J, T))
  for(i in 1:R) {
        M <- rpois(1, lambda) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
  }
  y <- matrix(y, nrow=R) # convert array to matrix

  ysc <- data.frame(cov1=rnorm(R*T))
  ysc$cov1[T] <- NA

  # Organize data
  umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
                          dist.breaks=breaks,
                          yearlySiteCovs=ysc,
                          tlength=rep(transect.length, R), numPrimary=T)

  # Fit the model
  expect_warning(fm1 <- gdistsamp(~1, ~1, ~cov1, umf, keyfun="exp", output="density", se=FALSE))

  # Check that getP works
  expect_warning(gp <- getP(fm1))
  expect_equivalent(dim(gp), c(R, T*J))
})

test_that("gdistsamp simulate method works",{

    set.seed(343)
    R <- 30
    T <- 3
    strip.width <- 50
    transect.length <- 200 #Area != 1
    breaks <- seq(0, 50, by=10)

    covs <- as.data.frame(matrix(rnorm(R*T),ncol=T))
    names(covs) <- paste0('par',1:3)

    beta <- c(0.4,0.3,0.6)
    lambda <- exp(1.3 + beta[1]*covs$par1)
    phi <- plogis(as.matrix(0.4 + beta[2]*covs))
    sigma <- exp(as.matrix(3 + beta[3]*covs))
    J <- length(breaks)-1
    y <- array(0, c(R, J, T))
    for(i in 1:R) {
        M <- rpois(1, lambda[i]) # Individuals within the 1-ha strip
        for(t in 1:T) {
            # Distances from point
            d <- runif(M, 0, strip.width)
            # Detection process
            if(length(d)) {
                cp <- phi[i,t]*exp(-d^2 / (2 * sigma[i,t]^2)) # half-normal w/ g(0)<1
                d <- d[rbinom(length(d), 1, cp) == 1]
                y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
    y <- matrix(y, nrow=R) # convert array to matrix

    covs$par1[2] <- NA
    umf <- unmarkedFrameGDS(y = y, siteCovs=covs, yearlySiteCovs=covs,
                            survey="line", unitsIn="m",
                            dist.breaks=breaks,
                            tlength=rep(transect.length, R), numPrimary=T)

    expect_warning(fm <- gdistsamp(~par1, ~1, ~1, umf, se=FALSE, K=15, engine="C"))

    #This used to error due to rmultinom not accepting size=NA
    expect_warning(s <- simulate(fm, nsim=2, na.rm=FALSE))
    expect_equivalent(length(s), 2)
    expect_equivalent(dim(s[[1]]), c(30,15))
    expect_true(!any(is.na(s[[1]][1,])))
    expect_true(all(is.na(s[[1]][2,])))

    expect_warning(pb <- parboot(fm, nsim=1))
    expect_is(pb, "parboot")

})

