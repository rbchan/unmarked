context("occu fitting function")

skip_on_cran()

test_that("occu can fit simple models",{

  y <- matrix(rep(1,10)[1:10],5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occu(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  expect_equivalent(occ,1)

  det <- coef(backTransform(det))
  expect_equivalent(det,1)

  bt <- backTransform(fm, type = 'state')
  expect_equivalent(coef(bt), 1)

  bt <- backTransform(fm, type = 'det')
  expect_equivalent(coef(bt), 1)

  est_obj <- fm@estimates@estimates$state
  expect_equal(est_obj@invlink, "logistic")
  expect_equal(est_obj@invlinkGrad, "logistic.grad")

  y <- matrix(rep(0,10)[1:10],5,2)
  umf <- unmarkedFrameOccu(y = y)
  fm <- occu(~ 1 ~ 1, data = umf)

  occ <- fm['state']
  det <- fm['det']

  occ <- coef(backTransform(occ))
  expect_equivalent(occ, 0, tolerance = 1e-4)

  det <- coef(backTransform(det))
  expect_equivalent(det,0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'state')
  expect_equivalent(coef(bt), 0, tolerance = 1e-4)

  bt <- backTransform(fm, type = 'det')
  expect_equivalent(coef(bt), 0, tolerance = 1e-4)


})

test_that("occu can fit models with covariates",{
  skip_on_cran()
  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ x, data = umf)
  fmR <- occu(~ o1 + o2 ~x, data = umf, engine="R")
  expect_equal(coef(fm), coef(fmR))

  occ <- fm['state']
  det <- fm['det']

  expect_error(occ <- coef(backTransform(occ)))

  expect_equivalent(coef(occ), c(8.590737, 2.472220), tolerance = 1e-4)
  expect_equivalent(coef(det), c(0.44457, -0.14706, 0.44103), tolerance = 1e-4)

  ci <- confint(occ)
  expect_equal(dim(ci), c(2,2))

  out <- capture.output(occ)
  expect_equal(out[1], "Occupancy:")

  occ.lc <- linearComb(fm, type = 'state', c(1, 0.5))
  det.lc <- linearComb(fm, type = 'det', c(1, 0.3, -0.3))

  expect_equivalent(coef(occ.lc), 9.826848, tol = 1e-4)
  expect_equivalent(coef(det.lc), 0.2681477, tol = 1e-4)

  expect_equivalent(coef(backTransform(occ.lc)), 1, tol = 1e-4)
  expect_equivalent(coef(backTransform(det.lc)), 0.5666381, tol = 1e-4)

  expect_error(backTransform(fm, type = "state"))
  expect_error(backTransform(fm, type = "det"))

  fitted <- fitted(fm)
  expect_equivalent(fitted, structure(c(0.5738, 0.5014, 0.4318, 0.38581, 0.50171, 0.53764,
0.46563, 0.40283, 0.39986, 0.79928), .Dim = c(5L, 2L)), tol = 1e-5)

  # methods
  gp <- getP(fm)
  expect_equal(dim(gp), c(5,2))
  res <- residuals(fm)
  expect_equal(dim(res), c(5,2))
  expect_equal(res[1,1], -0.57380, tol=1e-4)

  r <- ranef(fm)
  expect_equal(dim(r@post), c(5,2,1))
  expect_equal(bup(r), c(1,1,1,1,1))

  s <- simulate(fm, 2)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(umf@y))

  fitstats <- function(fm) {
      observed <- getY(fm@data)
      expected <- fitted(fm)
      resids <- residuals(fm)
      sse <- sum(resids^2,na.rm=TRUE)
      chisq <- sum((observed - expected)^2 / expected,na.rm=TRUE)
      freeTuke <- sum((sqrt(observed) - sqrt(expected))^2,na.rm=TRUE)
      out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
      return(out)
  }
  pb <- parboot(fm, fitstats, nsim=3)
  expect_equal(dim(pb@t.star), c(3,3))

  y <- matrix(rep(0,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  expect_warning(fm <- occu(~ o1 + o2 ~ x, data = umf))
  detMat <- fm@estimates@estimates$det@covMat
  stMat <- fm@estimates@estimates$state@covMat
  expect_equivalent(detMat, matrix(rep(NA,9),nrow=3))
  expect_equivalent(stMat, matrix(rep(NA,4),nrow=2))
  
  # Trap attempts to use a variable in formula that isn't in covariates
  fake <- rnorm(3)
  expect_error(fm <- occu(~ o1 + o2 ~ fake, data = umf))
  expect_error(fm <- occu(~ o1 + fake ~ o1, data = umf))
})

test_that("occu handles NAs",{

  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  siteCovs[3,1] <- NA
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  expect_warning(fm <- occu(~ o1 + o2 ~ x, data = umf))
  expect_equal(fm@sitesRemoved, 3)
  expect_equivalent(coef(fm), c(8.70123, 4.58255, 0.66243, -0.22862, 0.58192), tol = 1e-5)

  obsCovs[10,2] <- NA
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  expect_warning(fm <- occu(~ o1 + o2 ~ x, data = umf))
  expect_equal(fm@sitesRemoved, 3)
  expect_equivalent(coef(fm), c(8.91289, 1.89291, -1.42471, 0.67011, -8.44608), tol = 1e-5)

})

## Add some checks here.
test_that("occu handles offsets",{

  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ offset(x), data = umf)
  expect_equivalent(coef(fm),
                     structure(c(9.74361, 0.44327, -0.14683, 0.44085), .Names = c("psi(Int)",
"p(Int)", "p(o1)", "p(o2)")), tol = 1e-5)
  fm <- occu(~ o1 + offset(o2) ~ offset(x), data = umf)
  expect_equivalent(coef(fm), structure(c(8.59459, 0.97574, -0.3096), .Names = c("psi(Int)",
"p(Int)", "p(o1)")), tol=1e-5)

})

test_that("occu cloglog link function works",{
  skip_on_cran()
  #Adapted from example by J. Cohen
  set.seed(123)
  M = 500 #sample size
  J = 3 #number of visits

  #standardized covariates
  elev <- runif(n = M, 0,100)
  forest <- runif(n = M, 0,1)
  wind <- array(runif(n = M * J, 0,20), dim = c(M, J))
  elev=as.numeric(scale(elev))
  forest=as.numeric(scale(forest))
  wind[,1] <- as.numeric(scale(wind[,1]))
  wind[,2] <- as.numeric(scale(wind[,2]))
  wind[,3] <- as.numeric(scale(wind[,3]))

  #regression parameters for abundance
  beta0 = -0.69
  beta1 = 0.71
  beta2 = -0.5

  #simulate abundance and derive true occupancy
  lambda <- exp(beta0 + beta1 * elev + beta2 * forest)
  N <- rpois(n = M, lambda = lambda)
  z <- as.numeric(N>0)
  #regression parameters for detection
  alpha0 = -0.84
  alpha1 = 2.
  alpha2 = -1.2

  #simulate detection
  p <- plogis(alpha0 + alpha1 * elev + alpha2 * wind )

  #create vectors of simulation values, for easy comparison to model estimates
  true.beta.p <- c(alpha0,alpha1,alpha2)
  true.beta.occ <- c(beta0,beta1,beta2)

  #generate observed presence
  Obs.pres <- matrix(NA,M,J)
  for (i in 1:M){
    for (j in 1:J){
      Obs.pres[i,j] <- rbinom(1,1,z[i]*p[i,j])
    }
  }
  Obs.ever <- apply(Obs.pres,1,max)

  #create observation-level covariate data frame for unmarked
  sitevec <- rep(1:M,3) #vector of site ID's
  wind.df <- data.frame("wind"=wind)
  colnames(wind.df) <- c("Wind.1","Wind.2","Wind.3")
  wind.vec <- c(wind.df$Wind.1,wind.df$Wind.2,wind.df$Wind.3)
  wind.frame <- data.frame("site"=sitevec,"wind"=wind.vec)
  wind.frame.order <- wind.frame[order(wind.frame$site),]
  wind.for.um <- data.frame(wind.frame.order$wind)
  colnames(wind.for.um)="wind"

  #create unmarked data object
  occ.frame <- unmarkedFrameOccu(Obs.pres,
                                 siteCovs=data.frame("ele"=elev,"forest"=forest),
                                  obsCovs=wind.for.um)

  #create model object
  occ_test <-occu(~ele+wind ~ele+forest, occ.frame, linkPsi="cloglog",
                      se=F)
  truth <- c(true.beta.occ, true.beta.p)
  est <- coef(occ_test)
  expect_equivalent(truth, est, tol=0.1)
  expect_equivalent(est,
    c(-0.7425,0.6600,-0.3333,-0.87547,2.0677,-1.3082), tol=1e-4)

  est_obj <- occ_test@estimates@estimates$state
  expect_equal(est_obj@invlink, "cloglog")
  expect_equal(est_obj@invlinkGrad, "cloglog.grad")

  #Check error if wrong link function
  expect_error(occu(~ele+wind ~ele+forest, occ.frame, linkPsi="fake"))
})

test_that("occu predict works",{
  skip_on_cran()
  skip_if(!requireNamespace("raster", quietly=TRUE), 
          "raster package unavailable")
  set.seed(55)
  R <- 20
  J <- 4
  x1 <- rnorm(R)
  x2 <- factor(c(rep('A', R/2), rep('B', R/2)))
  x3 <- matrix(rnorm(R*J), R, J)
  z <- rbinom(R, 1, 0.5)
  y <- matrix(rbinom(R*J, 1, z*0.6), R, J)
  x1[1] <- NA
  x3[2,1] <- NA
  x3[3,] <- NA
  umf1 <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x1=x1, x2=x2),
                              obsCovs=list(x3=x3))
  fm1 <- expect_warning(occu(~x3 ~x1+x2, umf1))
  E1.1 <- expect_warning(predict(fm1, type="state"))
  E1.2 <- expect_warning(predict(fm1, type="det"))

  nd1.1 <- data.frame(x1=0, x2=factor('A', levels=c('A','B')))
  nd1.2 <- data.frame(x3=0)
  E1.3 <- predict(fm1, type="state", newdata=nd1.1, appendData=TRUE)
  E1.4 <- predict(fm1, type="det", newdata=nd1.2)

  r1 <- raster::raster(matrix(rnorm(100), 10))
  expect_error(predict(fm1, type="state", newdata=r1))
  s1 <- raster::stack(r1)
  expect_error(predict(fm1, type="state", newdata=s1))
  names(s1) <- c("x3")
  E1.5 <- predict(fm1, type="det", newdata=s1)
  E1.5 <- predict(fm1, type="det", newdata=s1, appendData=TRUE)

  E1.6 <- expect_warning(predict(fm1, type="state", level=0.9))
  expect_equal(as.numeric(E1.6[1,3:4]), c(0.01881844, 0.8538048))
})

test_that("occu predict can handle complex formulas",{

  y <- matrix(rep(0:1,10)[1:10],5,2)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10, o2 = exp(-5:4)/10)
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ scale(o1) + o2 ~ x, data = umf)

  #Predict values should not depend on means/variance of newdata itself
  nd1 <- obsCovs(umf[1:2,])
  pr1 <- predict(fm, 'det', newdata=nd1)
  nd2 <- obsCovs(umf[1:4,])
  pr2 <- predict(fm, 'det', newdata=nd2)[1:4,]

  expect_equivalent(pr1, pr2)

  #Check factors
  siteCovs$fac_cov <- factor(sample(c('a','b','c'), 5, replace=T),
                             levels=c('b','a','c'))

  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ fac_cov, data = umf)

  pr3 <- predict(fm, 'state', newdata=data.frame(fac_cov=c('a','b')))
  pr4 <- predict(fm, 'state', newdata=data.frame(fac_cov=c('b','a')))

  expect_equivalent(as.matrix(pr3),as.matrix(pr4[2:1,]))
  expect_error(predict(fm, 'state', newdata=data.frame(fac_cov=c('a','d'))))

  #Check when original covs contain factor not used in formula
  siteCovs$fac_cov2 <- factor(sample(c('a','b','c'), 5, replace=T),
                             levels=c('b','a','c'))

  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- occu(~ o1 + o2 ~ fac_cov, data = umf)
  #expect_warning(pr <- predict(fm, 'state', newdata=data.frame(fac_cov=c('a','b'))))

})


test_that("occu can handle random effects",{
  skip_on_cran()
  set.seed(123)
  n_sites <- 100
  n_years <- 8
  site_id <- rep(1:n_sites, each=n_years)
  M <- n_sites * n_years
  J <- 5 # number of obs per year
  site_covs <- data.frame(cov1=rnorm(M), site_id=factor(site_id))
  beta <- c(intercept=0.5, cov1=0.3)
  sig <- 1.2
  site_effect <- rnorm(n_sites, 0, sig)
  true_site_means <- plogis(beta[1] + site_effect)

  psi <- rep(NA, M)
  for (i in 1:M){
    #Random group intercept on psi
    psi[i] <- plogis(beta[1] + beta[2]*site_covs$cov1[i]
                      + site_effect[site_id[i]])
  }

  p <- 0.5
  z <- rbinom(M, 1, psi)
  y <- matrix(0, nrow=M, ncol=J)

  for (i in 1:M){
    if(z[i]==1){
      y[i,] <- rbinom(J, 1, p)
    }
  }

  umf <- unmarkedFrameOccu(y=y, siteCovs=site_covs)
  fm <- occu(~1~cov1 + (1|site_id), umf)
  expect_equivalent(coef(fm), c(0.65293, 0.39965, -0.02822), tol=1e-4)
  expect_equivalent(sigma(fm)$sigma, 1.18816, tol=1e-4)

  out <- capture.output(fm)
  expect_equal(out[6], "Random effects:")

  pr <- predict(fm, "state", newdata=data.frame(cov1=0, site_id=factor(1:100)))
  expect_is(pr, "data.frame")

  ft <- fitted(fm)
  expect_equivalent(dim(ft), c(n_sites*n_years, J))

  pb <- parboot(fm, nsim=1)
  expect_is(pb, "parboot")

  # confint should only show fixed effects
  ci <- confint(fm, type = 'state')
  expect_equal(nrow(ci), 2)

  ci <- confint(fm['state'])
  expect_equal(nrow(ci), 2)

  # Check custom initial values
  expect_equal(fm@TMB$starts_order[1], "beta_det")
  fmi <- occu(~1~cov1 + (1|site_id), umf, starts=c(10,0,0,0))
  expect_equivalent(fmi@TMB$par["beta_det"], 10)
  expect_error(occu(~1~cov1 + (1|site_id), umf, starts=rep(0,3)))
  expect_error(occu(~1~cov1 + (1|site_id), umf, starts=c(100,0,0,0)))

  # Check site covs as random effects in obs model
  fm <- occu(~(1|site_id)~1, umf)
  expect_true(sigma(fm)$Model[1]=="p")
  pr <- predict(fm, 'det')
  expect_true(inherits(pr, 'data.frame'))

  umf2 <- unmarkedFrameOccu(y=getY(umf), siteCovs=NULL,
  obsCovs=data.frame(obs_id=factor(sample(letters[1:5], length(getY(umf)), replace=T))))
  fm <- occu(~(1|obs_id)~1, umf2)
  expect_true(sigma(fm)$Model[1]=="p")

  # Check vcov method
  v1 <- vcov(fm)
  expect_equivalent(dim(v1), c(2,2))
  expect_equal(colnames(v1), c("psi(Int)","p(Int)"))

  v2 <- vcov(fm, fixedOnly=FALSE)
  expect_equivalent(dim(v2), c(7,7))
  expect_equal(colnames(v2), c("psi(Int)","p(Int)",rep("p(b_det)", 5)))

  fl <- fitList(m1=fm, m2=fm)
  #options(warn=2)
  #on.exit(options(warn=0))
  test <- modSel(fl) # shouldn't warn
  #options(warn=0)
})

test_that("TMB engine gives correct det estimates when there are lots of NAs", {

  skip_on_cran()
  set.seed(123)
  M <- 200
  J <- 10
  psi <- 0.5

  z <- rbinom(M, 1, psi)

  x <- matrix(rnorm(M*J), M, J)

  p <- plogis(0 + 0.3*x)

  y <- matrix(NA, M, J)
  for (i in 1:M){
    y[i,] <- rbinom(J, 1, p[i,]) * z[i]
  }
  y[sample(1:(M*J), (M*J/2), replace=FALSE)] <- NA

  umf <- unmarkedFrameOccu(y=y, obsCovs=list(x=x))

  fit <- occu(~x~1, umf)
  fitT <- occu(~x~1, umf, engine="TMB")
  expect_equal(coef(fit), coef(fitT), tol=1e-5)
})
