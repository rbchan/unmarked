test.pcount.offest <- function() 
{

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
  checkEqualsNumeric(coef(fm), structure(c(-0.78814924, 2.62569034, -0.02578801), 
      .Names = c("lam(Int)", "p(Int)", "p(o1)")), tol = 1e-5)

}
