context("ranef predict method")
skip_on_cran()

test_that("ranef predict method works",{
  #Single-season model
  set.seed(4564)
  R <- 10
  J <- 5
  N <- rpois(R, 3)
  y <- matrix(NA, R, J)
  y[] <- rbinom(R*J, N, 0.5)
  y[1,] <- NA
  y[2,1] <- NA
  K <- 15

  umf <- unmarkedFramePCount(y=y)
  fm <- expect_warning(pcount(~1 ~1, umf, K=K))

  re <- expect_warning(ranef(fm))

  set.seed(123)
  ps <- posteriorSamples(re, nsim=10)
  expect_is(ps, "unmarkedPostSamples")

  sh <- capture.output(show(ps))
  expect_equal(sh[1], "Posterior samples from unmarked model")

  #One is dropped bc of NA
  expect_equivalent(dim(ps@samples), c(9,1,10))

  # Brackets
  expect_equal(ps[1,1,1], ps@samples[1,1,1,drop=FALSE])

  # Method for unmarkedFit objects
  set.seed(123)
  ps2 <- expect_warning(posteriorSamples(fm, nsim=10))
  expect_equal(ps, ps2)

  # Custom function
  set.seed(123)
  myfunc <- function(x){
    c(gr1=mean(x[1:4]), gr2=mean(x[5:9]))
  }

  pr <- predict(re, fun=myfunc, nsim=10)
  expect_equivalent(dim(pr), c(2,10))
  expect_equal(rownames(pr), c("gr1","gr2"))
  expect_equivalent(as.numeric(pr[1,1:3]), c(7.0,5.0,4.75))

  #Dynamic model
  set.seed(7)
  M <- 10
  J <- 3
  T <- 5
  lambda <- 5
  gamma <- 0.4
  omega <- 0.6
  p <- 0.5
  N <- matrix(NA, M, T)
  y <- array(NA, c(M, J, T))
  S <- G <- matrix(NA, M, T-1)
  N[,1] <- rpois(M, lambda)
  y[,,1] <- rbinom(M*J, N[,1], p)
  for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        y[,,t+1] <- rbinom(M*J, N[,t+1], p)
  }

  # Prepare data
  umf <- unmarkedFramePCO(y = matrix(y, M), numPrimary=T)

  # Fit model and backtransform
  m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=20)
  re1 <- ranef(m1)

  ps <- posteriorSamples(re1, nsim=10)
  expect_equivalent(dim(ps@samples), c(10,5,10))
  expect_equivalent(ps@samples[1,,1],c(7,4,3,1,1))

  myfunc <- function(x){
    apply(x, 2, function(x) c(mean(x[1:4]), mean(x[5:9])))
  }

  pr <- predict(re1, fun=myfunc, nsim=10)
  expect_equivalent(dim(pr), c(2,5,10))
  expect_equivalent(pr[1,1:3,1], c(3.5,2.5,1.5))
})
