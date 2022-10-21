context("pcount.spHDS fitting function")
skip_on_cran()

test_that("pcount.spHDS can fit models",{
# Unit test here
set.seed(12)
x<- rnorm(50)
N<- rpois(50, exp(x) )
sigma<- 2
d<- seq(0.2, 4,,50)

p<- 1-exp(-exp( -(d^2)/(2*sigma*sigma) ) )
y<- rbinom(50, N, prob=p)

umf <- unmarkedFramePCount(y=matrix(y,ncol=1),
   siteCovs=data.frame(dist=d,Habitat=x))
#summary(umf)

fm1 <- pcount.spHDS(~ -1 + I(dist^2) ~ Habitat, umf, K=20)

expect_equivalent(
                   coef(fm1),
structure(c(-0.0521147712416163, 0.952296442491614, -1.66812493149504
), .Names = c("lam(Int)", "lam(Habitat)", "p(I(dist^2))"))
                    , tol=1e-5)

expect_error(pcount.spHDS(~(1|dummy)~1,umf))
})

