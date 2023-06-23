context("powerAnalysis method")
skip_on_cran()

test_that("powerAnalysis method works",{
  forms <- list(state=~elev, det=~1)
  coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
  design <- list(M=300, J=8) # 300 sites, 8 occasions per site
  occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)

  template_model <- occu(~1~elev, occu_umf)
  nul <- capture.output(expect_error(powerAnalysis(template_model)))

  effect_sizes <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

  nul <- capture_output({

  set.seed(123)
  pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05, nsim=10)
  expect_is(pa, "unmarkedPower")
  s <- summary(pa)$Power
  expect_true(s[2]>0.7)

  # output printout
  out <- capture.output(pa)
  expect_equal(out[5], "Power Statistics:")

  # update
  pa_up <- update(pa, alpha=0.5)
  expect_is(pa_up, "unmarkedPower")

  # fewer sites
  set.seed(123)
  pa2 <- powerAnalysis(template_model, effect_sizes, design=list(M=50, J=3), nsim=10)
  expect_true(all(summary(pa2)$Power < 0.4))

  # more sites
  set.seed(123)
  pa3 <- powerAnalysis(template_model, effect_sizes, design=list(M=400, J=4), nsim=10)
  expect_true(summary(pa3)$Power[2] >0.7)

  # set null
  set.seed(123)
  nul <- list(state=c(intercept=5, elev=0), det=c(intercept=0))
  pa4 <- powerAnalysis(template_model, effect_sizes, nulls=nul, nsim=10)
  expect_true(summary(pa4)$Power[1]==1)
  expect_equivalent(summary(pa4)$Null, c(5,0,0))

  # list
  pl <- unmarkedPowerList(list(pa, pa2, pa3, pa4))
  expect_is(pl, "unmarkedPowerList")
  s <- summary(pl)
  expect_is(s, "data.frame")

  pdf(NULL)
  pl_plot <- plot(pl)
  expect_is(pl_plot,"list")
  dev.off()

  # generate list
  scenarios <- expand.grid(M=c(50,100), J=c(2,3))
  pl <- unmarkedPowerList(template_model, effect_sizes, design=scenarios, nsim=10)
  expect_is(pl, "unmarkedPowerList")

  # With random effect
  set.seed(123)
  rguide <- list(group=factor(levels=letters[1:20]))
  rform <- list(state=~x+(1|group), det=~1)
  rcf <- list(state=c(intercept=0, x=0.5, group=0.7), det=c(intercept=0))
  umfr <- simulate("occu", formulas=rform, design=design, coefs=rcf, guide=rguide)
  fm <- occu(~1~x+(1|group), umfr)
  pa5 <- powerAnalysis(fm, rcf, nsim=10)
  s <- summary(pa5)
  expect_equal(nrow(s), 3)
  expect_equal(s$Power[2], 1)
  })
})

test_that("custom datasets can be passed to powerAnalysis",{
  set.seed(123)
  coefs <- list(state=c(intercept=0, elev=-0.3), det=c(intercept=0))
  design <- list(M=300, J=8)
  forms <- list(state=~elev, det=~1)
  pco_umf <- simulate("pcount", formulas=forms, coefs=coefs, design=design, nsim=10)

  # convert pcount to occu umf
  conv_umf <- lapply(pco_umf, function(x){
    y <- x@y
    y[y > 0] <- 1
    unmarkedFrameOccu(y=y, siteCovs=siteCovs(x),
                    obsCovs=obsCovs(x))
  })

  fit <- occu(~1~elev, conv_umf[[1]])

  nul <- capture.output({

  pa <- powerAnalysis(fit, coefs=coefs, datalist=conv_umf, nsim=10)
  expect_equivalent(summary(pa)$Power[2], 0.9, tol=1e-4)

  pa2 <- powerAnalysis(fit, coefs=coefs, nsim=10)
  expect_equivalent(summary(pa2)$Power[2], 0.8, tol=1e-4)

  })

  expect_error(powerAnalysis(fit, coefs=coefs, datalist=pco_umf))
  expect_error(powerAnalysis(fit, coefs=coefs, datalist=conv_umf, nsim=20))
})

test_that("powerAnalysis can be run in parallel",{
  skip_on_cran()
  skip_on_ci()
  forms <- list(state=~elev, det=~1)
  coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
  design <- list(M=50, J=3) # 300 sites, 8 occasions per site
  occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)

  template_model <- occu(~1~elev, occu_umf)
  nul <- capture.output(expect_error(powerAnalysis(template_model)))

  effect_sizes <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
  set.seed(123)
  pa <- powerAnalysis(template_model, coefs=effect_sizes, alpha=0.05, nsim=3,
                      parallel=TRUE)
  expect_is(pa, "unmarkedPower")


})
