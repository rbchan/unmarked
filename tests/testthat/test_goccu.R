context("goccu fitting function")
skip_on_cran()

set.seed(123)
M <- 100
T <- 5
J <- 4

psi <- 0.5
phi <- 0.3
p <- 0.4

z <- rbinom(M, 1, psi)
zmat <- matrix(z, nrow=M, ncol=T)

zz <- rbinom(M*T, 1, zmat*phi)
zz <- matrix(zz, nrow=M, ncol=T)

zzmat <- zz[,rep(1:T, each=J)]
y <- rbinom(M*T*J, 1, zzmat*p)
y <- matrix(y, M, J*T)
umf <- unmarkedMultFrame(y=y, numPrimary=T)

test_that("unmarkedFrameGOccu can be constructed", {
  set.seed(123)
  sc <- data.frame(x=rnorm(M))
  ysc <- matrix(rnorm(M*T), M, T)
  oc <- matrix(rnorm(M*T*J), M, T*J)

  umf2 <- unmarkedFrameGOccu(y, siteCovs=sc, obsCovs=list(x2=oc),
                           yearlySiteCovs=list(x3=ysc), numPrimary=T)
  expect_is(umf2, "unmarkedFrameGOccu")
  expect_equal(names(umf2@yearlySiteCovs), "x3")
})

test_that("goccu can fit models", {

  # Without covariates
  mod <- goccu(~1, ~1, ~1, umf)
  expect_equivalent(coef(mod), c(0.16129, -0.97041, -0.61784), tol=1e-5)
  
  # With covariates
  set.seed(123)
  sc <- data.frame(x=rnorm(M))
  ysc <- matrix(rnorm(M*T), M, T)
  oc <- matrix(rnorm(M*T*J), M, T*J)

  umf2 <- unmarkedMultFrame(y=y, siteCovs=sc, yearlySiteCovs=list(x2=ysc),
                            obsCovs=list(x3=oc), numPrimary=T)

  mod2 <- goccu(~x, ~x2, ~x3, umf2)
  expect_equivalent(coef(mod2), c(0.18895, -0.23629,-0.97246,-0.094335,-0.61808,
                                  -0.0040056), tol=1e-5)

  # predict
  pr <- predict(mod2, 'psi')
  expect_equal(dim(pr), c(M, 4))
  expect_equal(pr$Predicted[1], 0.5796617, tol=1e-5)

  # phi should not drop last level
  pr2 <- predict(mod2, 'phi')
  expect_equal(dim(pr2), c(M*T, 4))

  nd <- data.frame(x=1)
  pr3 <- predict(mod2, 'psi', newdata=nd)
  expect_true(nrow(pr3) == 1)
  expect_equal(pr3$Predicted[1], 0.488168, tol=1e-5)

  # Other methods
  ft <- fitted(mod2)
  expect_equal(dim(ft), dim(umf2@y))
  expect_true(all(ft >=0 & ft <= 1))

  res <- residuals(mod2)
  expect_equal(dim(res), dim(umf2@y))

  gp <- getP(mod2)
  expect_equal(dim(gp), dim(umf2@y))
  expect_equal(gp[1,1], 0.349239, tol=1e-5) 

  set.seed(123)
  s <- simulate(mod2, nsim=2)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(mod2@data@y))
  simumf <- umf2
  simumf@y <- s[[1]]
  simmod <- update(mod2, data=simumf)
  expect_equivalent(coef(simmod),
               c(0.174991, -0.27161, -1.32766, 0.054459,-0.41610,-0.073922), tol=1e-5)
  
  r <- ranef(mod2)
  expect_equal(dim(r@post), c(M, 2, 1))
  expect_equal(sum(bup(r)), 53.13565, tol=1e-4)
 
  pb <- parboot(mod2, nsim=2)
  expect_is(pb, "parboot")

  npb <- nonparboot(mod2, B=2, bsType='site')
  

})

test_that("goccu handles missing values", {

  set.seed(123)
  y2 <- y
  y2[1,1] <- NA
  y2[2,1:J] <- NA

  sc <- data.frame(x=rnorm(M))
  sc$x[3] <- NA
  ysc <- matrix(rnorm(M*T), M, T)
  ysc[4,1] <- NA
  oc <- matrix(rnorm(M*T*J), M, T*J)
  oc[5,1] <- NA
  oc[6,1:J] <- NA

  umf_na <- unmarkedMultFrame(y=y2, siteCovs=sc, yearlySiteCovs=list(x2=ysc),
                            obsCovs=list(x3=oc), numPrimary=T)
 
  mod_na <- expect_warning(goccu(~x, ~x2, ~x3, umf_na))
  
  pr <- expect_warning(predict(mod_na, 'psi'))
  expect_equal(nrow(pr), M-1)

  # Need to re-write these to use the design matrix instead of predict
  gp <- getP(mod_na)
  expect_equal(dim(gp), c(100, 20))
  expect_true(is.na(gp[5,1]))
  expect_true(all(is.na(gp[6, 1:4])))
  s <- simulate(mod_na)
  expect_equal(dim(s[[1]]), dim(mod_na@data@y))
  ft <- fitted(mod_na)
  expect_equal(dim(ft), dim(mod_na@data@y))
  r <- ranef(mod_na)
  expect_equal(dim(r@post), c(100, 2, 1))
  expect_true(is.na(bup(r)[3]))
 
  pb <- expect_warning(parboot(mod_na, nsim=2))
  expect_is(pb, "parboot")
})
