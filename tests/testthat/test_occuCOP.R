context("occuCOP fitting function")
skip_on_cran()

COPsimul <- function(psi = 0.5,
                     lambda = 1,
                     M = 100,
                     J = 5) {
  
  z_i <- sample(
    x = c(0, 1),
    size = M,
    prob = c(1 - psi, psi),
    replace = T
  )
  
  y = matrix(rpois(n = M * J, lambda = lambda), nrow = M, ncol = J) * z_i
  
  return(y)
}


test_that("unmarkedFrameOccuCOP is constructed correctly", {
  set.seed(123)
  
  # Simulate data
  M = 100
  J = 5
  y = COPsimul(psi = .5,
                   lambda = 1,
                   M = M,
                   J = J)
  L = y * 0 + 1
  
  psiCovs = data.frame(
    "psiNum" = rnorm(M),
    "psiInt" = as.integer(rpois(n = M, lambda = 5)),
    "psiBool" = sample(c(T, F), size = M, replace = T),
    "psiChar" = sample(letters[1:5], size = M, replace = T),
    "psiFactUnord" = factor(sample(
      letters[5:10], size = M, replace = T
    )),
    "psiFactOrd" = sample(factor(c("A", "B", "C"), ordered = T), size =
                            M, replace = T)
  )
  
  lambdaCovs = list(
    "lambdaNum" = matrix(
      rnorm(M * J), 
      nrow = M, ncol = J
    ),
    "lambdaInt" = matrix(
      as.integer(rpois(n = M * J, lambda = 1e5)),
      nrow = M, ncol = J
    ),
    "lambdaBool" = matrix(
      sample(c(T, F), size = M * J, replace = T),
      nrow = M, ncol = J
    ),
    "lambdaChar" = matrix(
      sample(letters[1:5], size = M * J, replace = T),
      nrow = M, ncol = J
    ),
    "lambdaFactUnord" = matrix(
      factor(sample(letters[5:10], size = M * J, replace = T)),
      nrow = M, ncol = J
    ),
    "lambdaFactOrd" = matrix(
      sample(factor(c("A", "B", "C"), ordered = T), size = M * J, replace = T),
      nrow = M, ncol = J
    )
  )

  
  # Creating a unmarkedFrameOccuCOP object
  expect_warning(umf <- unmarkedFrameOccuCOP(y = y))
  expect_no_error(umf <- unmarkedFrameOccuCOP(y = y, L = L))
  
  
  # Create subsets
  expect_no_error(umf_sub_i <- umf[1:3, ])
  expect_no_error(umf_sub_j <- umf[, 1:2])
  expect_no_error(umf_sub_ij <- umf[1:3, 1:2])
  
  # unmarkedFrameOccuCOP organisation ----------------------------------------------
  expect_true(inherits(umf, "unmarkedFrameOccuCOP"))
  expect_equivalent(numSites(umf_sub_i), 3)
  expect_equivalent(obsNum(umf_sub_j), 2)
  expect_equivalent(numSites(umf_sub_ij), 3)
  expect_equivalent(obsNum(umf_sub_ij), 2)
  
  # unmarkedFrameOccuCOP display ---------------------------------------------------
  
  # print method
  expect_output(print(umf), "Data frame representation of unmarkedFrame object")
  expect_output(print(umf_sub_i), "Data frame representation of unmarkedFrame object")
  expect_output(print(umf[1,]), "Data frame representation of unmarkedFrame object")
  expect_output(print(umf[,1]), "Data frame representation of unmarkedFrame object")
  expect_output(print(umf[1,1]), "Data frame representation of unmarkedFrame object")
  
  # summary method for unmarkedFrameOccuCOP
  expect_output(summary(umf), "unmarkedFrameOccuCOP Object")
  expect_output(summary(umf_sub_ij), "unmarkedFrameOccuCOP Object")
  
  # plot method for unmarkedFrameOccuCOP
  expect_no_error(plot(umf))
  expect_no_error(plot(umf_sub_ij))
  
  
  # Input handling: covariates -------------------------------------------------
  expect_no_error(umfCovs <- unmarkedFrameOccuCOP(
    y = y,
    L = L,
    siteCovs = psiCovs,
    obsCovs = lambdaCovs
  ))
  expect_output(print(umfCovs), "Data frame representation of unmarkedFrame object")
  expect_output(summary(umfCovs), "unmarkedFrameOccuCOP Object")
  expect_no_error(plot(umfCovs))
  
  # Input handling: NA ---------------------------------------------------------
  
  # NA should not pose problems when creating the unmarkedFrameOccuCOP object.
  # The warnings and potential errors will be displayed when fitting the model.
  # Except when y only contains NA: then there's an error.
  
  ## NA in y
  yNA <- y
  yNA[1:2,] <- NA
  yNA[3:5, 3:4] <- NA
  yNA[,3] <- NA
  expect_no_error(umfNA <- unmarkedFrameOccuCOP(y = yNA, L = L))
  expect_output(print(umfNA), "Data frame representation of unmarkedFrame object")
  expect_output(summary(umfNA), "unmarkedFrameOccuCOP Object")
  expect_no_error(plot(umfNA))
  
  ## NA in L
  obsLengthNA <- L
  obsLengthNA[1:2,] <- NA
  obsLengthNA[3:5, 3:4] <- NA
  obsLengthNA[,5] <- NA
  expect_no_error(umfNA <- unmarkedFrameOccuCOP(y = y, L = obsLengthNA))
  expect_output(print(umfNA), "Data frame representation of unmarkedFrame object")
  expect_output(summary(umfNA), "unmarkedFrameOccuCOP Object")

  expect_no_error(plot(umfNA))
  
  ## NA also in covariates
  psiCovsNA <- psiCovs
  psiCovsNA[1:5,] <- NA
  psiCovsNA[c(8,10,12), 3] <- NA
  psiCovsNA[,1] <- NA
  lambdaCovsNA <- lambdaCovs
  lambdaCovsNA[[1]][1:5,] <- NA
  lambdaCovsNA[[1]][,3] <- NA
  lambdaCovsNA[[3]][,5] <- NA
  expect_no_error(umfCovsNA <- unmarkedFrameOccuCOP(
    y = yNA,
    L = obsLengthNA,
    siteCovs = psiCovsNA,
    obsCovs = lambdaCovsNA
  ))
  expect_output(print(umfCovsNA), "Data frame representation of unmarkedFrame object")
  expect_output(summary(umfCovsNA), "unmarkedFrameOccuCOP Object")
  expect_no_error(plot(umfCovsNA))
  
  ## all NA in y
  yallNA <- y
  yallNA[1:M, 1:J] <- NA
  expect_error(unmarkedFrameOccuCOP(y = yallNA, L = L))
  
  # Input handling: error and warnings -----------------------------------------
  # Error when there are decimals in y
  y_with_decimals = y
  y_with_decimals[1, 1] = .5
  expect_error(unmarkedFrameOccuCOP(y = y_with_decimals, L = L))
  
  # Warning if y is detection/non-detection instead of count
  y_detec_nodetec = (y > 0) * 1
  expect_warning(unmarkedFrameOccuCOP(y = y_detec_nodetec, L = L))
  
  # Error if the dimension of L is different than that of y
  expect_error(unmarkedFrameOccuCOP(y = y, L = L[1:2, 1:2]))
})


test_that("occuCOP can fit simple models", {
  # Simulate data
  set.seed(123)
  M = 100
  J = 5
  y = COPsimul(psi = .5,
                   lambda = 1,
                   M = M,
                   J = J)
  L = y * 0 + 1

  # Create umf
  umf <- unmarkedFrameOccuCOP(y = y, L = L)

  # Fitting options ----

  ## With default parameters ----
  expect_no_error(fit_default <- occuCOP(data = umf, L1 = TRUE))
  expect_warning(occuCOP(data = umf, psiformula = ~ 1, lambdaformula = ~ 1, psistarts = 0, lambdastarts = 0))

  ## With chosen starting points ----
  expect_no_error(occuCOP(data = umf,
                          psiformula = ~ 1, lambdaformula = ~ 1,
                          psistarts = qlogis(.7),
                          lambdastarts = log(0.1), L1=T))
  expect_no_error(occuCOP(data = umf,
                          psiformula = ~ 1, lambdaformula = ~ 1,
                          starts = c(qlogis(.7), log(0.1)), L1 = T))
  # warning if all starts and psistarts and lambdastarts were furnished
  # and starts != c(psistarts, lambdastarts)
  expect_no_error(occuCOP(data = umf, starts = c(0, 0),
                          psistarts = c(0), lambdastarts = c(0), L1 = T))
  expect_warning(occuCOP(data = umf, starts = c(0, 1),
                          psistarts = c(0), lambdastarts = c(0), L1 = T))
  # errors if starting vectors of the wrong length
  expect_error(occuCOP(data = umf, starts = c(0), L1 = T))
  expect_error(occuCOP(data = umf, psistarts = c(0, 0), lambdastarts = 0, L1 = T))
  expect_error(occuCOP(data = umf, lambdastarts = c(0, 0), L1 = T))

  # With different options
  expect_no_error(occuCOP(data = umf, method = "Nelder-Mead", psistarts = 0, lambdastarts = 0, L1=T))
  expect_error(occuCOP(data = umf, method = "ABC", psistarts = 0, lambdastarts = 0, L1=T))

  expect_no_error(occuCOP(data = umf, se = F, psistarts = 0, lambdastarts = 0, L1=T))
  expect_error(occuCOP(data = umf, se = "ABC"))

  expect_no_error(occuCOP(data = umf, engine = "R", psistarts = 0, lambdastarts = 0, L1=T))
  expect_error(occuCOP(data = umf, engine = "julia", psistarts = 0, lambdastarts = 0, L1=T))

  expect_no_error(occuCOP(data = umf, na.rm = F, psistarts = 0, lambdastarts = 0, L1=T))
  expect_error(occuCOP(data = umf, na.rm = "no", psistarts = 0, lambdastarts = 0, L1=T))

  # Looking at at COP model outputs ----
  expect_is(fit_default, "unmarkedFitOccuCOP")
  expect_equivalent(coef(fit_default), c(0.13067954, 0.06077929), tol = 1e-5)
  
  ## backTransform
  expect_no_error(backTransform(fit_default, type = "psi"))
  expect_no_error(backTransform(fit_default, type = "lambda"))
  expect_error(backTransform(fit_default, type = "state"))
  expect_error(backTransform(fit_default, type = "det"))
  expect_is(backTransform(fit_default, type = "psi"), "unmarkedBackTrans")
  
  ## predict with newdata = fit@data
  expect_no_error(umpredpsi <- predict(object = fit_default, type = "psi"))
  expect_equal(umpredpsi$Predicted[1], 0.5326235, tol = 1e-5)
  expect_no_error(umpredlambda <- predict(object = fit_default, type = "lambda"))
  expect_equal(umpredlambda$Predicted[1], 1.062664, tol = 1e-5)
  expect_error(predict(object = fit_default, type = "state"))
  
  ## predict with newdata = 1
  expect_no_error(predict(
    object = fit_default,
    newdata = data.frame(1),
    type = "psi"
  ))
  expect_no_error(predict(
    object = fit_default,
    newdata = data.frame(1),
    type = "lambda"
  ))
  expect_no_error(predict(
    object = fit_default,
    newdata = data.frame("a"=1:5,"b"=10:14),
    type = "psi"
  ))
  
  # Fitting accurately ----
  ## psi = 0.50, lambda = 1 ----
  psi_test = .5
  lambda_test = 1
  fit_accur <- occuCOP(data = unmarkedFrameOccuCOP(
    y = COPsimul(
      psi = psi_test,
      lambda = lambda_test,
      M = 1000,
      J = 10
    ),
    L = matrix(1, nrow = 1000, ncol = 10)
  ), psistarts = 0, lambdastarts = 0, L1=T)
  psi_estimate = backTransform(fit_accur, type = "psi")@estimate
  lambda_estimate = backTransform(fit_accur, type = "lambda")@estimate
  expect_equivalent(
    psi_estimate,
    psi_test,
    tol = 0.05
  )
  expect_equivalent(
    lambda_estimate,
    lambda_test,
    tol = 0.05
  )

  ## psi = 0.25, lambda = 5 ----
  psi_test = 0.25
  lambda_test = 5
  fit_accur <- occuCOP(data = unmarkedFrameOccuCOP(
    y = COPsimul(
      psi = psi_test,
      lambda = lambda_test,
      M = 1000,
      J = 10
    ),
    L = matrix(1, nrow = 1000, ncol = 10)
  ), psistarts = 0, lambdastarts = 0, L1=T)
  psi_estimate = backTransform(fit_accur, type = "psi")@estimate
  lambda_estimate = backTransform(fit_accur, type = "lambda")@estimate
  expect_equivalent(
    psi_estimate,
    psi_test,
    tol = 0.05
  )
  expect_equivalent(
    lambda_estimate,
    lambda_test,
    tol = 0.05
  )

  ## psi = 0.75, lambda = 0.5 ----
  psi_test = 0.75
  lambda_test = 0.5
  fit_accur <- occuCOP(data = unmarkedFrameOccuCOP(
    y = COPsimul(
      psi = psi_test,
      lambda = lambda_test,
      M = 1000,
      J = 10
    ),
    L = matrix(1, nrow = 1000, ncol = 10)
  ), psistarts = 0, lambdastarts = 0, L1=T)
  psi_estimate = backTransform(fit_accur, type = "psi")@estimate
  lambda_estimate = backTransform(fit_accur, type = "lambda")@estimate
  expect_equivalent(
    psi_estimate,
    psi_test,
    tol = 0.05
  )
  expect_equivalent(
    lambda_estimate,
    lambda_test,
    tol = 0.05
  )

  # With NAs ----
  yNA <- y
  yNA[1,] <- NA
  yNA[3, 1] <- NA
  yNA[4, 3] <- NA
  yNA[, 5] <- NA
  expect_no_error(umfNA <- unmarkedFrameOccuCOP(y = yNA, L = L))

  expect_warning(fit_NA <- occuCOP(data = umfNA, psistarts = 0, lambdastarts = 0, L1=T))
  expect_error(occuCOP(data = umfNA, psistarts = 0, lambdastarts = 0, na.rm = F))
})

test_that("We can simulate COP data", {
  
  # From scratch ----
  
  # With no covariates
  expect_no_error(simulate(
    "occuCOP",
    formulas = list(psi =  ~ 1, lambda =  ~ 1),
    coefs = list(
      psi = c(intercept = 0),
      lambda = c(intercept = 0)
    ),
    design = list(M = 100, J = 100)
  ))
  
  # With quantitative covariates
  expect_no_error(simulate(
    "occuCOP",
    formulas = list(psi =  ~ elev, lambda =  ~ rain),
    coefs = list(
      psi = c(intercept = qlogis(.5), elev = -0.5),
      lambda = c(intercept = log(3), rain = -1)
    ),
    design = list(M = 100, J = 5)
  ))
  
  # With guides
  expect_no_error(simulate(
    "occuCOP",
    formulas = list(psi =  ~ elev, lambda =  ~ rain),
    coefs = list(
      psi = c(intercept = qlogis(.5), elev = -0.5),
      lambda = c(intercept = log(3), rain = -1)
    ),
    design = list(M = 100, J = 5),
    guide = list(elev=list(dist=rnorm, mean=12, sd=0.5))
  ))
  
  # With qualitative covariates
  expect_no_error(umf <- simulate(
    "occuCOP",
    formulas = list(psi =  ~ elev + habitat, lambda =  ~ 1),
    coefs = list(
      psi = c(
        intercept = qlogis(.2),
        elev = -0.5,
        habitatB = .5,
        habitatC = .8
      ),
      lambda = c(intercept = log(3))
    ),
    design = list(M = 100, J = 5),
    guide = list(habitat = factor(levels = c("A", "B", "C")))
  ))
  
  # From unmarkedFitOccuCOP ----
  expect_no_error(umfit <- occuCOP(
    umf,
    psiformula =  ~ habitat,
    L1 = T,
    psistarts = c(0,0,0),
    lambdastarts = 0
  ))
  expect_no_error(simulate(umfit))
})

test_that("occuCOP can fit and predict models with covariates", {
  # Simulate data with covariates ----
  set.seed(123)
  expect_no_error(umf <- simulate(
    "occuCOP",
    formulas = list(psi =  ~ elev + habitat, lambda =  ~ rain),
    coefs = list(
      psi = c(
        intercept = qlogis(.2),
        elev = -0.5,
        habitatB = .5,
        habitatC = .8
      ),
      lambda = c(intercept = log(3), rain = -1)
    ),
    design = list(M = 100, J = 5),
    guide = list(habitat = factor(levels = c("A", "B", "C")))
  ))
  
  # Fit ----
  expect_no_error(umfit <- occuCOP(
    umf,
    psiformula =  ~ habitat + elev,
    lambdaformula =  ~ rain,
    L1 = T,
    psistarts = c(0,0,0,0),
    lambdastarts = c(0,0)
  ))
  
  expect_error(occuCOP(
    umf,
    psiformula =  ~ habitat+elev,
    lambdaformula =  ~ rain,
    L1 = T,
    psistarts = c(0),
    lambdastarts = c(0,0)
  ))
  
  expect_equivalent(
    coef(umfit), 
    c(-1.5350679, 0.4229763, 0.7398768, -1.0456397, 1.2333424, -0.8344109), 
    tol = 1e-5
  )
  
  # Predict ----
  expect_no_error(predict(umfit, type = "psi"))
  expect_no_error(umpredpsi <- predict(
    umfit,
    type = "psi",
    newdata = data.frame("habitat" = c("A", "B", "C"), "elev" = c(0, 0, 0)),
    appendData = TRUE
  ))
  expect_equivalent(umpredpsi$Predicted, c(0.1772534, 0.2474811, 0.3110551), tol = 1e-5)
  
  expect_no_error(umpredlambda <- predict(umfit, type = "lambda", appendData = TRUE))
  expect_no_error(predict(umfit, type = "lambda", level = 0.5))
  expect_equal(umpredlambda$Predicted[1], 1.092008, tol = 1e-5)
})

