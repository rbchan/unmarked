context("occuFP fitting function")

test_that("occuFP model can be fit",{
  n = 100
  o = 10
  o1 = 5
  y = matrix(0,n,o)
  p = .7
  r = .5
  fp = 0.05
  y[1:(n*.5),(o-o1+1):o] <- rbinom((n*o1*.5),1,p)
  y[1:(n*.5),1:(o-o1)] <- rbinom((o-o1)*n*.5,1,r)
  y[(n*.5+1):n,(o-o1+1):o] <- rbinom((n*o1*.5),1,fp)
  type <- c((o-o1),o1,0)
  site <- c(rep(1,n*.5*.8),rep(0,n*.5*.2),rep(1,n*.5*.2),rep(0,n*.8*.5))
  occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
  site <- data.frame(habitat = site)
  occ <- list(METH = occ)
  umf1 <- unmarkedFrameOccuFP(y,site,occ, type = type)

  m1 <- occuFP(detformula = ~ METH, FPformula = ~1,
               stateformula = ~ habitat, data = umf1)
  expect_equal(names(m1), c("state","det","fp"))
  expect_warning(fl <- fitList(m1,m1))
  expect_is(fl,"unmarkedFitList")
  expect_equal(length(fl@fits), 2)

  pr <- predict(m1, "fp")
  expect_equal(dim(pr), c(1000, 4))

  # Check error when random effect in formula
  expect_error(occuFP(~(1|dummy), ~1, ~1, data=umf1))
})
