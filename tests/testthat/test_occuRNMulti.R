context("occuRNMulti fitting function")
skip_on_cran()

# Simulate data
set.seed(1)
M <- 300
J <- 5
N1 <- N2 <- N3 <- 25
while((max(N1) > 20) | (max(N2) > 20) | (max(N3) > 20)){
# Covariates
sc <- data.frame(x1 = rnorm(M), x2 = rnorm(M))
oc <- data.frame(x3 = rnorm(M*J))
# Coyote expected abundance
l1 <- 2
# True coyote abundance
N1 <- rpois(M, l1)
# Lambda for marten
l2 <- exp(log(1.3) - 0.3*N1 + 0.5*sc$x1*N1)
# True marten abundance
N2 <- rpois(M, l2)
# Lambda for fisher
l3 <- exp(log(1) + -0.2 * N2)
N3 <- rpois(M, l3)
}
# True detection prob for individuals of both species
rs <- 0.2
ps1 <- 1 - (1-rs)^N1
ps2 <- 1 - (1-rs)^N2
ps3 <- 1 - (1-rs)^N3
# True parameter values
truth <- c(log(2), log(1.3), -0.3, 0.5, log(1), -0.2, rep(log(0.2/(1-0.2)), 3)) 
# Simulate observations
y1 <- y2 <- y3 <- matrix(NA, M, J)
for (m in 1:M){
  for (j in 1:J){
    y1[m,] <- rbinom(J, 1, ps1[m])
    y2[m,] <- rbinom(J, 1, ps2[m])
    y3[m,] <- rbinom(J, 1, ps3[m])
  }
}
umf <- unmarkedFrameOccuMulti(y=list(coyote=y1, marten=y2, fisher=y3),
                              siteCovs=sc, obsCovs=oc)

test_that("occuRNMulti can fit 2-species models",{
  umf2 <- umf
  umf2@ylist <- umf@ylist[1:2]
  sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~x1))
  df <- list(coyote=~1, marten=~1)
  fit <- occuRNMulti(df, sf, umf2[1:50,])
  expect_is(fit, "unmarkedFitOccuRNMulti") 
  expect_equivalent(coef(fit), c(0.3888, 0.4018,-0.4229,0.6316,-1.0916,-1.8663),
                    tol=1e-4)
  
  pr <- suppressMessages(predict(fit, type='state'))
  expect_equal(length(pr), 2)
  expect_equal(nrow(pr[[1]]), 50)

  nd <- data.frame(x1=c(0,1))
  pr <- suppressMessages(predict(fit, type='state', newdata=nd))
  expect_equal(pr[[2]]$Predicted, c(0.8009, 2.0334), tol=1e-4)

  res <- residuals(fit)
  expect_equal(length(res), 2)
  expect_equal(res[[1]][1], -0.30977, tol=1e-4)

  s <- simulate(fit, nsim=2)
  expect_equal(length(s), 2)
  expect_equal(length(s[[1]]), 2)

  pb <- parboot(fit, nsim=2)
  expect_is(pb, "parboot")

  # Unsupported methods
  expect_error(ranef(fit))
  expect_error(nonparboot(fit))
  expect_error(linearComb(fit, type='state', c(0,1,1)))

  # Check with NAs
  umf2@ylist[[1]][1,1] <- NA
  umf2@ylist[[2]][1,1] <- NA
  umf2@ylist[[1]][2,] <- NA
  umf2@ylist[[2]][2,] <- NA
  umf2@siteCovs$x1[3] <- NA # missing covariate value

  fit_na <- occuRNMulti(df, sf, umf2[1:50,])
  expect_equal(coef(fit), coef(fit_na), tol=0.2)

  pr <- suppressMessages(predict(fit_na, type='state'))
  expect_true(is.na(pr[[2]]$Predicted[3]))

  df2 <- list(coyote=~x1, marten=~1)
  fit_na2 <- occuRNMulti(df2, sf, umf2[1:20,])
  pr <- predict(fit_na2, type='det')
  expect_true(all(is.na(pr$coyote$Predicted[11:15])))

  # Occupancy model for subordinate
  fit_occ <- occuRNMulti(df, sf, umf2, modelOccupancy='marten')
  expect_equal(names(fit_occ), c("lam", "psi", "det"))
  pr <- suppressMessages(predict(fit_occ, type='state'))
  expect_true(max(pr$marten$Predicted, na.rm=TRUE) < 1)
})

test_that("occuRNMulti can fit 3-species models", {
            
  # `sp1 --> sp2 --> sp3`
  sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~x1),
           fisher = list(~1, marten = ~1))
  df <- list(coyote=~1, marten=~1, fisher=~1)

  fit <- occuRNMulti(df, sf, umf[1:10,])
  expect_equal(length(coef(fit)), 9)
  pr <- suppressMessages(predict(fit, 'state'))
  expect_equal(length(pr), 3)
  
  # `sp2 <-- sp1 --> sp3`
  sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~x1),
           fisher = list(~1, coyote = ~1))
  
  fit <- occuRNMulti(df, sf, umf[1:10,])
  expect_equal(length(coef(fit)), 9)
  expect_equal(
  names(coef(fit, type='lam')), 
  c("lam([coyote] (Intercept))", "lam([fisher] (Intercept))", 
    "lam([fisher:coyote] (Intercept))",
    "lam([marten] (Intercept))", "lam([marten:coyote] (Intercept))",
    "lam([marten:coyote] x1)")
  )
  pr <- suppressMessages(predict(fit, 'state'))
  expect_equal(length(pr), 3)

  # `sp1 --> sp3 <-- sp2`
  sf <- list(coyote = ~1,
           marten = ~1,
           fisher = list(~1, coyote = ~1, marten = ~1))
  fit <- occuRNMulti(df, sf, umf[1:10,])
  expect_equal(names(coef(fit, type='lam')),
    c("lam([coyote] (Intercept))", "lam([marten] (Intercept))", 
      "lam([fisher] (Intercept))", "lam([fisher:coyote] (Intercept))", 
      "lam([fisher:marten] (Intercept))"
      ))
  pr <- suppressMessages(predict(fit, 'state'))
  expect_equal(length(pr), 3)

  # `sp1 --> sp2 --> sp3`
  # plus sp1 --> sp3
  sf <- list(coyote = ~1,
           marten = list(~1, coyote = ~x1),
           fisher = list(~1, coyote=~1, marten = ~1))
  df <- list(coyote=~1, marten=~1, fisher=~1)

  # This model fit doesn't work well due to tiny sample size
  fit <- suppressWarnings(occuRNMulti(df, sf, umf[1:10,]))
  expect_equal(names(coef(fit, 'lam')),
  c("lam([coyote] (Intercept))", "lam([marten] (Intercept))", 
    "lam([marten:coyote] (Intercept))", "lam([marten:coyote] x1)", 
    "lam([fisher] (Intercept))", "lam([fisher:coyote] (Intercept))",
    "lam([fisher:marten] (Intercept))"))
  pr <- suppressMessages(predict(fit, 'state', level=NULL))
  expect_equal(length(pr), 3)

})
