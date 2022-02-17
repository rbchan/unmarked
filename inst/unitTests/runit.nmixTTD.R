#Setup common vars
M = 100 # Number of sites
nrep <- 3 # Number of visits per site
Tmax = 5 # Max duration of a visit
alpha1 = -1 # Covariate on rate
beta1 = 1 # Covariate on density
mu.lambda = 1 # Rate at alpha1 = 0
mu.dens = 1 # Density at beta1 = 0

set.seed(123)
covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
covDens <- rnorm(M) #Abundance/density covariate

test.nmixTTD.P.exp <- function(){
  set.seed(123)
  dens <- exp(log(mu.dens) + beta1 * covDens)
  N <- rpois(M, dens) # Realized density per site
  lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
  ttd <- NULL
  for(i in 1:nrep){
    ttd <- cbind(ttd,rexp(M, N*lambda[,i]))
  }
  ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
  ttd[ttd >= Tmax] <- 5 # Crop at Tmax
  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens,
                                                  cdens2=rnorm(length(covDens)),
                                                  cdens3=rnorm(length(covDens))),
                            obsCovs = data.frame(covDet=as.vector(t(covDet)),
                                                 cdet2=rnorm(length(covDet)),
                                                 cdet3=rnorm(length(covDet))))
  fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10)
  checkEqualsNumeric(coef(fit), c(-0.2846,1.1224,0.2221,-1.06713), tol=1e-4)
  #with NA
  umf@y[1,1] <- NA
  fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10)
  checkEqualsNumeric(coef(fit), c(-0.2846,1.1224,0.2221,-1.06713), tol=1e-4)

  #Predict
  pr1 <- predict(fit, "state")
  checkEqualsNumeric(dim(pr1), c(M, 4))
  checkEqualsNumeric(pr1$Predicted[1:2], c(0.3371,0.3232), tol=1e-4)
  nd <- data.frame(covDens=0)
  pr2 <- predict(fit, "state", newdata=nd)
  checkEqualsNumeric(pr2$Predicted, 0.7523, tol=1e-4)
  pr3 <- predict(fit, "det")
  checkEqualsNumeric(dim(pr3), c(M*nrep, 4))
  checkEqualsNumeric(pr3$Predicted[1], 2.2710, tol=1e-4)
  nd <- data.frame(covDet=0)
  pr4 <- predict(fit, "det", newdata=nd)
  checkEqualsNumeric(dim(pr4), c(1,4))
  checkEqualsNumeric(pr4$Predicted, 1.248748, tol=1e-4)

  #Check with site covs in det formula
  fit_pred <- nmixTTD(~1, ~covDens, data=umf, K=max(N)+10)
  pr5 <- predict(fit_pred, "det")
  checkEqualsNumeric(dim(pr5), c(M*nrep, 4))
  checkEqualsNumeric(pr5$Predicted[1:3], rep(0.587956, 3), tol=1e-4)
  checkEqualsNumeric(pr5$Predicted[4:6], rep(0.5769142, 3), tol=1e-4)
  nd5 <- data.frame(covDens=c(1,2))
  pr6 <- predict(fit_pred, "det", newdata=nd5)
  checkEqualsNumeric(dim(pr6), c(2,4))

  #Simulate
  sim <- simulate(fit, 2)
  checkTrue(inherits(sim, "list"))
  checkEqualsNumeric(length(sim), 2)
  checkEqualsNumeric(dim(sim[[1]]), dim(umf@y))

  #Update
  fit2 <- update(fit, data=umf[1:10,])

  #Ranef
  r <- ranef(fit)
  b <- bup(r)
  checkEqualsNumeric(length(b), M)
  checkEqualsNumeric(b[2], 1.204738, tol=1e-5)

  checkException(residuals(fit))

  #Try with threads=2
  fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10)
  fit_2 <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10, threads=2)
  checkEqualsNumeric(coef(fit), coef(fit_2))

  #Check error when random effect in formula
  checkException(nmixTTD(~(1|dummy), ~1, umf))

  #Check with more than 1 detection covariate
  fit3 <- nmixTTD(~covDens, ~covDet+cdet2+cdet3, data=umf, K=max(N)+10)
  checkEqualsNumeric(length(coef(fit3, "det")), 4)
  fit4 <- nmixTTD(~covDens+cdens2, ~covDet+cdet2+cdet3, data=umf, K=max(N)+10)
  checkEqualsNumeric(length(coef(fit4, "state")), 3)
  checkEqualsNumeric(length(coef(fit4, "det")), 4)
}

test.nmixTTD.P.weib <- function(){
  set.seed(123)
  shape = 5
  dens <- exp(log(mu.dens) + beta1 * covDens)
  N <- rpois(M, dens)
  lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
  ttd <- NULL
  for(i in 1:nrep) {
    ttd <- cbind(ttd,rweibull(M, shape, 1/(N*lambda[,i])))
  }
  ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
  ttd[ttd >= Tmax] <- 5 # Crop at Tmax
  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens),
                            obsCovs = data.frame(covDet=as.vector(t(covDet))))

  fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10, ttdDist="weibull")
  checkEquals(names(fit@estimates@estimates), c("state","det","shape"))
  checkEqualsNumeric(coef(fit), c(-0.08528,1.0540,0.0326,-0.9981,1.7203), tol=1e-4)

  sim <- simulate(fit, 2)
  r <- ranef(fit)
}

test.nmixTTD.NB.exp <- function(){
  set.seed(123)
  M = 500 # Number of sites
  covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
  covDens <- rnorm(M) #Abundance/density covariate
  dens <- exp(log(mu.dens) + beta1 * covDens)
  disp=2
  N <- rnbinom(M, mu=dens, size=disp) # Realized density per site
  lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
  ttd <- NULL
  for(i in 1:nrep) {
    ttd <- cbind(ttd,rexp(M, N*lambda[,i]))  # Simulate time to first detection per visit
  }
  ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
  ttd[ttd >= Tmax] <- 5 # Crop at Tmax
  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens),
                            obsCovs = data.frame(covDet=as.vector(t(covDet))))

  fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10, mixture="NB")
  checkEquals(names(fit@estimates@estimates), c("state","det","alpha"))
  checkEqualsNumeric(coef(fit), c(-0.0991,1.017,-0.00264,-0.99387, 1.1499), tol=1e-4)

  sim <- simulate(fit, 2)
  r <- ranef(fit)
}

test.nmixTTD.NB.weib <- function(){
  set.seed(123)
  M = 500 # Number of sites
  covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
  covDens <- rnorm(M) #Abundance/density covariate
  dens <- exp(log(mu.dens) + beta1 * covDens)
  disp=2
  shape = 5
  N <- rnbinom(M, mu=dens, size=disp) # Realized density per site
  lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
  ttd <- NULL
  for(i in 1:nrep) {
    ttd <- cbind(ttd,rweibull(M, shape, 1/(N*lambda[,i])))
  }
  ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
  ttd[ttd >= Tmax] <- 5 # Crop at Tmax
  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens),
                            obsCovs = data.frame(covDet=as.vector(t(covDet))))
  fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10,
                 mixture="NB", ttdDist="weibull")
  checkEquals(names(fit@estimates@estimates), c("state","det","alpha","shape"))
  checkEqualsNumeric(coef(fit), c(-0.0605,1.0328,-0.0006,-1.008,0.8141,1.6311), tol=1e-4)

  sim <- simulate(fit, 2)
  r <- ranef(fit)
}

test.nmixTTD.ZIP <- function(){

  M = 500 # Number of sites
  set.seed(123)
  covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
  covDens <- rnorm(M) #Abundance/density covariate

  dens <- exp(log(mu.dens) + beta1 * covDens)
  psi <- 0.3
  N <- unmarked:::rzip(M, dens, psi) # Realized density per site
  lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
  ttd <- NULL
  for(i in 1:nrep){
    ttd <- cbind(ttd,rexp(M, N*lambda[,i]))
  }
  ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
  ttd[ttd >= Tmax] <- 5 # Crop at Tmax
  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens,
                                                  cdens2=rnorm(length(covDens)),
                                                  cdens3=rnorm(length(covDens))),
                            obsCovs = data.frame(covDet=as.vector(t(covDet)),
                                                 cdet2=rnorm(length(covDet)),
                                                 cdet3=rnorm(length(covDet))))
  fit <- nmixTTD(~covDens, ~covDet, data=umf, mixture="ZIP", K=max(N)+10)
  checkEqualsNumeric(coef(fit), c(-0.0586,0.9033,0.1039,-1.0365,-0.9251), tol=1e-4)
  
  # Predict
  pr <- predict(fit, 'state')
  checkTrue(nrow(pr) == length(dens))
  checkException(predict(fit, 'psi'))
  bt <- backTransform(fit, 'psi')
  checkEquals(bt@estimate, 0.2839, tol=1e-4)

  # Simulate
  s <- simulate(fit, 2)
  checkTrue(length(s)== 2) 

  # ranef
  r <- ranef(fit)
  b <- bup(r)
  checkEquals(cor(b, N), 0.9148, tol=1e-4)

  # Compare R and C engines
  M = 50 # Number of sites
  set.seed(123)
  covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
  covDens <- rnorm(M) #Abundance/density covariate

  dens <- exp(log(mu.dens) + beta1 * covDens)
  psi <- 0.3
  N <- unmarked:::rzip(M, dens, psi) # Realized density per site
  lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
  ttd <- NULL
  for(i in 1:nrep){
    ttd <- cbind(ttd,rexp(M, N*lambda[,i]))
  }
  ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
  ttd[ttd >= Tmax] <- 5 # Crop at Tmax
  umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens,
                                                  cdens2=rnorm(length(covDens)),
                                                  cdens3=rnorm(length(covDens))),
                            obsCovs = data.frame(covDet=as.vector(t(covDet)),
                                                 cdet2=rnorm(length(covDet)),
                                                 cdet3=rnorm(length(covDet))))

  fit <- nmixTTD(~covDens, ~covDet, data=umf, mixture="ZIP", K=max(N)+10)
  fitR <- nmixTTD(~covDens, ~covDet, data=umf, mixture="ZIP", K=max(N)+10, engine="R")
  checkEqualsNumeric(coef(fit), coef(fitR))
}
