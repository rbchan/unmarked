context("IDS fitting function")
skip_on_cran()

test_that("IDS can fit models with covariates", {
  set.seed(123)
  formulas <- list(lam=~elev, ds=~1, phi=~1)
  # Based on values from real data
  design <- list(Mds=1000, J=6, Mpc=300)
  # Based on values from real data
  coefs <- list(lam = c(intercept=3, elev=-0.5),
              ds = c(intercept=-2.5),
              phi = c(intercept=-1.3))
  # Survey durations, loosely based on real data
  durs <- list(ds = rep(5, design$Mds), pc=runif(design$Mpc, 3, 30))

  sim_umf <- simulate("IDS", # name of model we are simulating for
                    nsim=1, # number of replicates
                    formulas=formulas, 
                    coefs=coefs,
                    design=design,
                    # arguments used by unmarkedFrameDS
                    dist.breaks = seq(0, 0.30, length.out=7),
                    unitsIn="km", 
                    # arguments used by IDS
                    # could also have e.g. keyfun here
                    durationDS=durs$ds, durationPC=durs$pc, durationOC=durs$oc,
                    maxDistPC=0.5, maxDistOC=0.5,
                    unitsOut="kmsq") 
  set.seed(123)
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1,
                dataDS=sim_umf$ds, dataPC=sim_umf$pc,
                availformula = ~1, durationDS=durs$ds, durationPC=durs$pc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")

  expect_equivalent(coef(mod_sim), c(3.0271179,-0.4858101,-2.5050244,-1.3729505), tol=1e-5)

  # Predict
  pr <- predict(mod_sim, type = 'lam')
  expect_equal(nrow(pr), 1000) # predicts only for the distance sampled sites
  pr <- predict(mod_sim, type = 'lam', newdata=sim_umf$pc)

  pr <- predict(mod_sim, type = 'ds')
  expect_equal(nrow(pr), 1000)
  expect_equal(pr$Predicted[1], exp(-2.5), tol=0.05)
  pr <- predict(mod_sim, type = 'pc')
  expect_equal(nrow(pr), 300)
  expect_equal(pr$Predicted[1], exp(-2.5), tol=0.05)
  pr <- predict(mod_sim, type = 'phi')
  expect_equal(pr$ds$Predicted[1], exp(-1.37), tol=0.05)
  
  # residuals
  r <- residuals(mod_sim)
  expect_equal(lapply(r, dim), list(ds=c(1000,6), pc = c(300,1)))

  # parboot
  pb <- parboot(mod_sim, nsim=2)
  
  pdf(NULL)
  plot(mod_sim)
  hist(mod_sim)
  dev.off()

  expect_error(nonparboot(mod_sim))
  expect_error(ranef(mod_sim))

  # Separate detection models
  mod_sep <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaPC = ~1,
                dataDS=sim_umf$ds[1:100,], dataPC=sim_umf$pc[1:100,],
                availformula = ~1, durationDS=durs$ds[1:100], durationPC=durs$pc[1:100],
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")
  expect_equal(length(coef(mod_sim)), 4)
  expect_equal(length(coef(mod_sep)), 5)

  # Doesn't support random effects
  expect_error(update(mod_sep, lambdaformula = ~elev + (1|group)))
})

test_that("IDS can fit models with occupancy data", {

  set.seed(123)
  formulas <- list(lam=~elev, ds=~1, oc=~1)
  # Based on values from real data
  design <- list(Mds=100, J=6, Mpc=100, Moc=100)
  # Based on values from real data
  coefs <- list(lam = c(intercept=3, elev=-0.5),
              ds = c(intercept=-2.5),
              oc = c(intercept = -2))

  sim_umf <- simulate("IDS", # name of model we are simulating for
                    nsim=1, # number of replicates
                    formulas=formulas, 
                    coefs=coefs,
                    design=design,
                    # arguments used by unmarkedFrameDS
                    dist.breaks = seq(0, 0.30, length.out=7),
                    unitsIn="km", 
                    # arguments used by IDS
                    # could also have e.g. keyfun here
                    maxDistPC=0.5, maxDistOC=0.5,
                    unitsOut="kmsq")

  mod_oc <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=sim_umf$ds, dataPC=sim_umf$pc, dataOC=sim_umf$oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")

  expect_equivalent(coef(mod_oc), c(3.0557091, -0.4200719, -2.5384331, -2.0610341),
                    tol=1e-5)
  
  pr <- predict(mod_oc, type='oc')
  expect_equal(pr$Predicted[1], 0.1273222, tol=1e-5)
  
  res <- residuals(mod_oc)
  expect_equal(lapply(res, dim), list(ds=c(100,6), pc = c(100,1), oc=c(100,1)))

  pb <- parboot(mod_oc, nsim=1)

  # Don't estimate availability if OC data
  expect_error(
  mod_oc <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                availformula = ~1,
                dataDS=sim_umf$ds, dataPC=sim_umf$pc, dataOC=sim_umf$oc,
                maxDistPC=0.5, maxDistOC=0.5,
                durationDS=durs$ds, durationPC=durs$pc, durationOC=durs$pc,
                unitsOut="kmsq")
  )

  # Just occupancy data
  mod_oc <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=sim_umf$ds, dataOC=sim_umf$oc,
                maxDistOC=0.5,
                unitsOut="kmsq")
  
  expect_equal(names(mod_oc), c("lam","ds","oc"))

})

test_that("IDS handles missing values", {

  set.seed(123)
  design <- list(Mds=100, J=6, Mpc=100, Moc=100)
  formulas <- list(lam=~elev, ds=~1, phi=~1)
  # Based on values from real data
  coefs <- list(lam = c(intercept=3, elev=-0.5),
              ds = c(intercept=-2.5),
              phi = c(intercept=-1.3))
  # Survey durations, loosely based on real data
  durs <- list(ds = rep(5, design$Mds), pc=runif(design$Mpc, 3, 30))

  sim_umf <- simulate("IDS", # name of model we are simulating for
                    nsim=1, # number of replicates
                    formulas=formulas, 
                    coefs=coefs,
                    design=design,
                    # arguments used by unmarkedFrameDS
                    dist.breaks = seq(0, 0.30, length.out=7),
                    unitsIn="km", 
                    # arguments used by IDS
                    # could also have e.g. keyfun here
                    maxDistPC=0.5, maxDistOC=0.5,
                    unitsOut="kmsq")

  sim_umf$pc@y[1,1] <- NA
  sim_umf$pc@y[2,] <- NA

  sim_umf$oc@y[1,1] <- NA
  sim_umf$oc@y[2,] <- NA

  expect_warning(
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=sim_umf$ds, dataPC=sim_umf$pc, dataOC=sim_umf$oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")
  )

  expect_equivalent(coef(mod_sim), c(2.9354934,-0.4759405,-2.5314594,-2.3259133),
                    tol=1e-5)

  
  sim_umf$ds@y[1,1] <- NA
  sim_umf$ds@y[2,] <- NA

  expect_error(
  expect_warning(
  mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1, detformulaOC = ~1,
                dataDS=sim_umf$ds, dataPC=sim_umf$pc, dataOC=sim_umf$oc,
                maxDistPC=0.5, maxDistOC=0.5,
                unitsOut="kmsq")
  ))

})
