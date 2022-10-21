context("pcount fitting function")
skip_on_cran()

test_that("pcount can fit simple models",{

  y <- matrix(c(
      8,7,
      6,7,
      8,8,
      8,6,
      7,7), nrow=5, ncol=2, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:10)
  umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- pcount(~ o1 ~ offset(x), data = umf, K=30)
  expect_equivalent(coef(fm), structure(c(-0.78814924, 2.62569034, -0.02578801),
      .Names = c("lam(Int)", "p(Int)", "p(o1)")), tol = 1e-5)

  y <- matrix(c(
      8,7,7,8,
      6,7,7,5,
      8,8,7,8,
      4,5,5,5,
      4,4,3,3), nrow=5, ncol=4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = seq(-1, 1, length=length(y)))
  umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fm <- pcount(~ o1 ~ x, data = umf, K=30)
  expect_equivalent(coef(fm),
                     c(1.91984184, -0.02987393,  2.49421875, -0.23350448),
                     tol = 1e-5)

  fm <- pcount(~o1~x, data=umf, K=30, mixture='ZIP')
  pr <- predict(fm, 'state')
  expect_is(pr, "data.frame")
  expect_equivalent(as.numeric(pr[1,]), c(6.8184,2.0456,3.7871,12.2760), tol=1e-4)

  nd <- data.frame(x=c(0,1))
  pr <- predict(fm, 'state', newdata=nd)
  expect_equal(dim(pr), c(2,4))

  pr <- predict(fm, 'det')
  expect_equal(dim(pr), c(20,4))

  gp <- getP(fm)
  expect_equal(dim(gp), dim(umf@y))

  res <- residuals(fm)
  expect_equal(dim(res), dim(umf@y))

  r <- ranef(fm)
  expect_equal(dim(r@post), c(5,31,1))
  expect_equal(bup(r), c(8.01, 7.01, 8.07, 5.03, 4.01), tol=1e-3)
  fm2 <- update(fm, mixture="NB")
  r2 <- ranef(fm2)
  expect_is(r2, "unmarkedRanef")
  fm3 <- update(fm, mixture="ZIP")
  r3 <- ranef(fm3)
  expect_is(r3, "unmarkedRanef")

  s <- simulate(fm, n=2)
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]]), dim(umf@y))

  pb <- parboot(fm, nsim=1)
  expect_is(pb, "parboot")

})

test_that("pcount predict works",{

  set.seed(55)
  R <- 20
  J <- 4
  N <- rpois(R, 2)
  y <- matrix(rbinom(R*J, N, 0.7), R, J)
  umf1 <- unmarkedFramePCount(y=y)

  fm1 <- pcount(~1 ~1, umf1, K=40)
  E1.1 <- predict(fm1, type="state")
  E1.2 <- predict(fm1, type="det")

  fm2 <- pcount(~1 ~1, umf1, K=40, mixture="NB")
  E2.1 <- predict(fm2, type="state")
  expect_error(predict(fm2, type="alpha"))

  fm3 <- pcount(~1 ~1, umf1, K=40, mixture="ZIP")
  E3.1 <- predict(fm3, type="state")
  expect_error(predict(fm3, type="psi"))
  expect_equal(E3.1[1,1], 1.818512, tol=1e-6)

})

test_that("pcount can fit models with random effects",{

set.seed(35)
nSites <- 300
nVisits <- 3
x <- rnorm(nSites)               # a covariate
beta0 <- 0
beta1 <- 0.4

ran <- rnorm(100, 0, 1)
group <- factor(as.character(rep(1:100, each=3)))
ran_ind <- as.numeric(group)

lambda <- exp(beta0 + beta1*x +
              ran[ran_ind])   # expected counts at each site
N <- rpois(nSites, lambda)       # latent abundance
y <- matrix(NA, nSites, nVisits)
p <- c(0.3, 0.6, 0.8)            # detection prob for each visit
for(j in 1:nVisits) {
  y[,j] <- rbinom(nSites, N, p[j])
}

# Organize data
visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)

expect_warning(umf <- unmarkedFramePCount(y=y, siteCovs=data.frame(x=x,group=group),
         obsCovs=list(visit=visitMat)))

fm <- pcount(~1~x, umf, K=25)
expect_is(fm, "unmarkedFitPCount")

fmr <- pcount(~visit~x+(1|group), umf, K=25)

expect_equivalent(coef(fmr), c(0.05397,0.3197,-0.8760,1.3668,2.078),
                   tol=1e-3)

expect_true(inherits(sigma(fmr), 'data.frame'))
expect_equal(sigma(fmr)$sigma, 1.05945, tol=1e-5)

pr <- predict(fmr, "state")
expect_equivalent(as.numeric(pr[1,]),
                   c(1.037050,0.58179,0.3453,3.1140), tol=1e-3)

pr2 <- predict(fmr, "state", re.form=NA)
expect_equivalent(as.numeric(pr2[1,]),
                   c(1.48366,0.2011,1.1374,1.93255), tol=1e-3)

pr3 <- predict(fmr, "det")
expect_true(inherits(pr3, "data.frame"))

nd <- data.frame(x=siteCovs(umf)$x[c(1,4)], group=factor(c(1,2)))
pr4 <- predict(fmr, "state", newdata=nd)
expect_equivalent(pr4$Predicted, pr$Predicted[c(1,4)])

# New group level
nd <- data.frame(x=c(0,1), group=factor(101))
expect_error(predict(fmr, "state", newdata=nd))

nd <- data.frame(x=c(0,1))
expect_error(expect_warning(predict(fmr, "state", newdata=nd)))

pr5 <- predict(fmr, "state", newdata=nd, re.form=NA)
expect_true(inherits(pr5, "data.frame"))

ft <- fitted(fmr)
expect_equivalent(dim(ft), c(300,3))

r <- ranef(fmr)
expect_true(inherits(r, "unmarkedRanef"))
b <- bup(r)
expect_true(cor(N, b) > 0.95)

rt <- randomTerms(fmr)
expect_true(inherits(rt, "data.frame"))
expect_equivalent(dim(rt), c(100,8))
expect_true(cor(ran, rt$Estimate) > 0.8)


# Multiple random effects
umf2 <- umf
siteCovs(umf2)$id <- sample(letters[1:3], 300, replace=T)

fmr2 <- pcount(~1~x+(1|group)+(1|id), umf2, K=25)

expect_true(nrow(sigma(fmr2))==2)
rt2 <- randomTerms(fmr2)
expect_true(all(rt2$Groups==c(rep("group",100), rep("id",3))))

# Check other distributions
fmnb <- pcount(~1~1, umf, engine="TMB", mixture="NB", K=25)
expect_true(inherits(fmnb@TMB, "list"))
expect_true(all(names(fmnb@estimates@estimates)==c("state","det","alpha")))

fmzip <- pcount(~1~1, umf, engine="TMB", mixture="ZIP", K=25)
expect_true(inherits(fmnb@TMB, "list"))
expect_true(all(names(fmnb@estimates@estimates)==c("state","det","alpha")))

# Site random effects in det formula
fm <- pcount(~(1|group)~1, umf2, K=19)
expect_true(sigma(fm)$Model[1]=="p")
})

test_that("pcount R, C++ and TMB engines give same results",{

  y <- matrix(c(
      8,7,7,8,
      6,7,7,5,
      8,8,7,8,
      4,5,5,5,
      4,4,3,3), nrow=5, ncol=4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = seq(-1, 1, length=length(y)))
  umf <- unmarkedFramePCount(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
  fmC <- pcount(~ o1 ~ x, data = umf, K=30, control=list(maxit=1))
  fmT <- pcount(~ o1 ~ x, data = umf, K=30, control=list(maxit=1), engine="TMB")
  fmR <- pcount(~ o1 ~ x, data = umf, K=30, control=list(maxit=1), engine="R")
  expect_equal(coef(fmC), coef(fmR))
  expect_equal(coef(fmC), coef(fmT), tol=1e-7)
})
