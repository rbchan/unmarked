context("occuMulti fitting function")
skip_on_cran()

test_that("unmarkedFrameOccuMulti construction and methods work",{

  y <- list(matrix(1:15,5,3),
            matrix(1:15,5,3))
  umf <- unmarkedFrameOccuMulti(y = y)

  expect_is(umf, "unmarkedFrameOccuMulti")
  out <- capture.output(umf)
  expect_equal(out[2], "Only showing observation matrix for species 1.")
  s <- capture.output(summary(umf))
  expect_equal(s[4], "2 species: sp1 sp2 ")

  # Check plot
  pdf(NULL)
  pl <- plot(umf)
  expect_is(pl, "trellis")
  dev.off()

  # check subset
  umf_sub <- umf[,1:2]
  expect_equal(umf_sub@ylist[[1]], umf@ylist[[1]][,1:2])
})

test_that("occuMulti can fit simple models",{

  y <- list(matrix(rep(1,10)[1:10],5,2),
            matrix(rep(1,10)[1:10],5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)

  #Probably should not be calling predict here b/c unit test
  #but complicated to get actual occupancy prob otherwise
  occ <- predict(fm,'state')$Predicted[1,1]
  expect_equivalent(occ,1, tolerance = 1e-4)

  detlist <- predict(fm,'det')
  det <- sapply(detlist,function(x) x[1,1])
  expect_equivalent(det, rep(1,length(detlist)), tolerance= 1e-4)

  #Check fitList
  expect_message(fl <- fitList(fm, fm))
  expect_message(expect_warning(fl <- fitList(fm, fm, autoNames='formula')))
  expect_message(expect_warning(fl <- fitList(fits=list(fm, fm), autoNames='formula')))
  expect_is(fl,"unmarkedFitList")
  expect_equivalent(length(fl@fits), 2)

  #Check error when random effect in formula
  expect_error(occuMulti(detformulas=rep("~1",2),
                           stateformulas=c("~(1|group)",rep("~1",2)), umf))

  y <- list(matrix(rep(0,10)[1:10],5,2),
            matrix(rep(0,10)[1:10],5,2))
  umf <- unmarkedFrameOccuMulti(y = y)
  fm <- occuMulti(detformulas=rep("~1",2),
                  stateformulas=rep("~1",3), data = umf, se=FALSE)

  occ <- predict(fm,'state')$Predicted[1,1]
  expect_equivalent(occ,0, tolerance = 1e-4)

  detlist <- predict(fm,'det')
  det <- sapply(detlist,function(x) x[1,1])
  expect_equivalent(det, rep(0,length(detlist)), tolerance= 1e-4)


})

test_that("occuMulti can fit models with covariates",{

  y <- list(matrix(rep(0:1,10)[1:10],5,2),
            matrix(rep(0:1,10)[1:10],5,2))

  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)
  stateformulas <- c('~occ_cov1','~occ_cov2','~occ_cov3')
  detformulas <- c('~det_cov1','~det_cov2')

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)

  occ <- fm['state']
  det <- fm['det']

  expect_equivalent(coef(occ), c(5.36630,0.79876,5.45492,-0.868451,9.21242,1.14561),
                     tolerance = 1e-4)
  expect_equivalent(coef(det), c(-0.27586,-0.81837,-0.09537,0.42334), tolerance = 1e-4)

  fit <- fitted(fm)
  expect_equivalent(length(fit),2)
  expect_equivalent(sapply(fit,function(x) x[1,1]),c(0.14954,0.30801), tol = 1e-4)

  res <- residuals(fm)
  expect_equivalent(length(res),2)
  expect_equivalent(sapply(res,function(x) x[1,1]),c(-0.14954,-0.30801), tol= 1e-4)

  gp <- getP(fm)
  expect_equivalent(length(gp), 2)
  expect_equivalent(dim(gp[[1]]), c(N,J))

  # ranef
  expect_error(ran <- ranef(fm))
  ran <- ranef(fm, species=1)
  expect_equal(bup(ran), rep(1,5))

  #Check site cov can be used in detection formula
  detformulas <- c('~occ_cov1','~det_cov2')
  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)
  expect_equivalent(coef(fm,'det')[2],3.355328e-05, tol=1e-4)
})

test_that("occuMulti can handle NAs",{

  y <- list(matrix(rep(0:1,10)[1:10],5,2),
            matrix(rep(0:1,10)[1:10],5,2))

  set.seed(456)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~occ_cov1','~occ_cov2','~occ_cov3')
  detformulas <- c('~det_cov1','~det_cov2')

  #Check error thrown when missing site covariates
  occ_covsNA <- occ_covs
  occ_covsNA[1,1] <- NA
  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covsNA, obsCovs = det_covs)
  expect_error(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))

  yna <- y
  yna[[1]][1,1] <- NA
  expect_warning(umf <- unmarkedFrameOccuMulti(y = yna, siteCovs = occ_covs, obsCovs = det_covs))

  #Check correct answer given when missing detection
  expect_warning(fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE))
  expect_equivalent(coef(fm)[c(1,7)], c(6.63207,0.35323), tol= 1e-4)

  fit <- fitted(fm)
  expect_true(is.na(fit[[1]][1,1]))

  res <- residuals(fm)
  expect_true(is.na(res[[1]][1,1]))

  gp <- getP(fm)
  expect_true(is.na(gp[[1]][1,1]))

  #Check error thrown when all detections are missing
  yna[[1]][1,] <- NA
  expect_warning(umf <- unmarkedFrameOccuMulti(y = yna, siteCovs = occ_covs, obsCovs = det_covs))
  expect_error(occuMulti(detformulas, stateformulas, data=umf, se=FALSE))

  #Check warning when missing covariate value on detection
  det_covsNA <- det_covs
  det_covsNA[1,1] <- NA
  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covsNA)
  expect_warning(occuMulti(detformulas,stateformulas,data=umf, se=FALSE))

})

test_that("occuMulti handles fixed 0 parameters",{

  y <- list(matrix(rep(0:1,10)[1:10],5,2),
            matrix(rep(0:1,10)[1:10],5,2))

  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~occ_cov1','~occ_cov2','0')
  detformulas <- c('~det_cov1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE)

  occ <- fm['state']
  expect_equivalent(length(coef(occ)),4)
  expect_equivalent(coef(occ),c(12.26043,0.61183,12.41110,0.18764),tol=1e-4)


  stateformulas <- c('~occ_cov1','~occ_cov2')
  fm2 <- occuMulti(detformulas, stateformulas, data = umf, maxOrder=1,se=FALSE)

  occ <- fm2['state']
  expect_equivalent(length(coef(occ)),4)
  expect_equivalent(coef(occ),c(12.26043,0.61183,12.41110,0.18764),tol=1e-4)

})

test_that("occuMulti predict method works",{

  set.seed(123)
  y <- list(matrix(rbinom(40,1,0.2),20,2),
            matrix(rbinom(40,1,0.3),20,2))

  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~occ_cov1','~occ_cov2','0')
  detformulas <- c('~det_cov1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  nul <- capture.output(prState <- predict(fm, type='state'))
  expect_equivalent(sapply(prState,function(x) x[1,1]),
                     c(0.30807707,0.20007250,0.04234835,0.73106618),tol=1e-4)
  nul <- capture.output(prDet <- predict(fm, type='det'))
  expect_equivalent(as.numeric(prDet$sp2[1,]),
                     c(0.190485,0.12201,0.0475270,0.525988), tol=1e-4)

  #Check with newdata
  nd <- siteCovs(umf)[1:2,]
  nul <- capture.output(pr_nd <- predict(fm, type='state', newdata=nd)$Predicted)
  expect_equivalent(pr_nd[,1],c(0.3080771,0.3196486), tol=1e-4)
  nd <- siteCovs(umf)[1:2,]
  nul <- capture.output(pr_nd <- predict(fm, type='state', newdata=nd, species=1, cond=2)$Predicted)
  expect_equivalent(pr_nd,c(0.3858233,0.5402935), tol=1e-4)
  #Make sure it works with newdata having only one row
  nd <- siteCovs(umf)[1,]
  nul <- capture.output(pr_nd <- predict(fm, type='state', newdata=nd)$Predicted)
  expect_equivalent(pr_nd[,1],c(0.3080771), tol=1e-4)
  nul <- capture.output(pr_nd <- predict(fm, type='state', newdata=nd, species=1, cond=2)$Predicted)
  expect_equivalent(pr_nd,c(0.3858233), tol=1e-4)

  stateformulas <- c('~1','~1','0')
  detformulas <- c('~1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  nul <- capture.output(prState <- predict(fm, type='state'))
  expect_equivalent(sapply(prState,function(x) x[1,1]),
                     c(0.475928,0.2548407,0.01496681,0.86713789),tol=1e-4)
  nul <- capture.output(prDet <- predict(fm, type='det'))
  expect_equivalent(as.numeric(prDet$sp2[1,]),
                     c(0.20494,0.11865,0.0582563,0.517888), tol=1e-4)

  #Check predicting co-occurrence

  nul <- capture_output({

  nd <- siteCovs(umf)[1:2,]
  pr_all <- predict(fm, type='state', se=F)$Predicted[1:2,1]
  pr_nd <- predict(fm, type='state', newdata=nd, species=c(1,2))$Predicted
  expect_equivalent(pr_nd,pr_all, tol=1e-4)

  #Check with site cov in detection formula
  stateformulas <- c('~occ_cov2','~1','0')
  detformulas <- c('~occ_cov1','~det_cov2')
  fm <- occuMulti(detformulas, stateformulas, data = umf)
  pr_state_actual <- predict(fm, "state")
  expect_equivalent(length(pr_state_actual), 4)
  expect_equivalent(pr_state_actual$Predicted[1,1], 0.729927907, tol=1e-5)
  expect_equivalent(nrow(pr_state_actual$Predicted), 20)

  pr_det_actual <- predict(fm, "det")
  expect_equivalent(length(pr_det_actual), 2)
  expect_equivalent(pr_det_actual$sp1$Predicted[1], 0.1448311, tol=1e-5)
  expect_equivalent(nrow(pr_det_actual$sp1), 20*2)

  #with newdata
  pr_state_nd <- predict(fm, "state", newdata=data.frame(occ_cov2=0))
  expect_equivalent(length(pr_state_nd), 4)
  expect_equivalent(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  expect_equivalent(nrow(pr_state_nd$Predicted), 1)

  pr_det_nd <- predict(fm, "det", newdata=data.frame(occ_cov1=0, det_cov2=0))
  expect_equivalent(length(pr_det_nd), 2)
  expect_equivalent(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  expect_equivalent(nrow(pr_state_nd$Predicted), 1)

  #With maxOrder set
  stateformulas <- c('~occ_cov2','~1')
  detformulas <- c('~occ_cov1','~det_cov2')

  fm <- occuMulti(detformulas, stateformulas, data = umf, maxOrder=1)

  pr_state <- predict(fm, "state")
  expect_equivalent(length(pr_state), 4)
  expect_equivalent(pr_state$Predicted[1,1], 0.729927907, tol=1e-5)
  expect_equivalent(nrow(pr_state$Predicted), 20)

  pr_state_nd <- predict(fm, "state", newdata=data.frame(occ_cov2=0))
  expect_equivalent(length(pr_state_nd), 4)
  expect_equivalent(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  expect_equivalent(nrow(pr_state_nd$Predicted), 1)

  pr_det <- predict(fm, "det")
  expect_equivalent(length(pr_det), 2)
  expect_equivalent(pr_det$sp1$Predicted[1], 0.1448311, tol=1e-5)
  expect_equivalent(nrow(pr_det$sp1), 20*2)

  pr_det_nd <- predict(fm, "det", newdata=data.frame(occ_cov1=0, det_cov2=0))
  expect_equivalent(length(pr_det_nd), 2)
  expect_equivalent(pr_state_nd$Predicted[1,1], 0.7538309, tol=1e-5)
  expect_equivalent(nrow(pr_state_nd$Predicted), 1)

  })

  #getP with maxOrder set
  gp <- getP(fm)
  expect_equal(length(gp), 2)
  expect_equal(dim(gp[[1]]), c(20,2))

  #simulate with maxOrder set
  s <- simulate(fm, 2)
  expect_is(s, "list")
  expect_equal(length(s), 2)
  expect_equal(dim(s[[1]][[1]]), c(N, J))

  #fitList with maxOrder set
  fm2 <- occuMulti(c("~1","~1"), c("~1","~1"), umf, maxOrder=1)
  expect_message(fl2 <- fitList(fm, fm2))
  expect_is(fl2, "unmarkedFitList")
  ms <- modSel(fl2)
  expect_is(ms, "unmarkedModSel")

  #fitted with maxOrder set
  ft <- fitted(fm)
  expect_equal(length(ft), 2)

  #parboot with maxOrder set
  pb <- parboot(fm, nsim=2)
  expect_is(pb, "parboot")


})

test_that("occuMulti predict can handle NAs",{

  set.seed(123)
  y <- list(matrix(rbinom(40,1,0.2),20,2),
            matrix(rbinom(40,1,0.3),20,2))

  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')
  det_covs[1,1] <- NA

  stateformulas <- c('~occ_cov1','~occ_cov2','0')
  detformulas <- c('~det_cov1','~det_cov2')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  expect_warning(fm <- occuMulti(detformulas, stateformulas, data = umf))

  prDet <- predict(fm, type='det')
  expect_true(all(is.na(prDet$sp1[1,])))
  expect_equivalent(as.numeric(prDet$sp1[2,]),
                     c(0.49781,0.19621,0.175514,0.8219401), tol=1e-4)

  #Check that you can predict with NAs in siteCovs

  nul <- capture_output({

  newdata <- siteCovs(umf)
  newdata[1,1] <- NA
  prOcc <- predict(fm, type='state', newdata=newdata)
  expect_true(all(is.na(prOcc$Predicted[1,])))
  expect_true(all(!is.na(sapply(prOcc,`[`,2,1))))
  prOcc_sp <- predict(fm, type='state', species=1, newdata=newdata)
  expect_true(all(is.na(prOcc_sp[1,])))
  expect_true(all(!is.na(prOcc_sp[2,])))
  expect_equivalent(prOcc_sp$Predicted[2],0.4731427, tol=1e-4)
  prOcc_cond <- predict(fm, type='state', species=1, cond=2, newdata=newdata)
  expect_true(all(is.na(prOcc_cond[1,])))
  expect_true(all(!is.na(prOcc_cond[2,])))
  expect_equivalent(prOcc_sp$Predicted[2],0.4731427, tol=1e-4)

  })
})


test_that("occuMulti can handle complex formulas",{

  #Check scale(), etc
  set.seed(123)
  y <- list(matrix(rbinom(40,1,0.2),20,2),
            matrix(rbinom(40,1,0.3),20,2))

  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3, mean=2),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2, mean=3),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  stateformulas <- c('~scale(occ_cov1)','~1','0')
  detformulas <- c('~scale(det_cov1)','~1')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)

  fm <- occuMulti(detformulas, stateformulas, data = umf)

  nul <- capture_output({

  #Check with newdata; contents of newdata should not
  #effect resulting predictions (scale should be based on
  #original data)
  nd <- siteCovs(umf)[1:5,]
  pr_nd <- predict(fm, type='state', newdata=nd, se=F)$Predicted
  nd <- siteCovs(umf)[1:2,]
  pr_nd2 <- predict(fm, type='state', newdata=nd, se=F)$Predicted
  nd <- siteCovs(umf)[c(1,1),]
  pr_nd3 <- predict(fm, type='state', newdata=nd, se=F)$Predicted

  expect_equivalent(pr_nd[1:2,], pr_nd2)
  expect_equivalent(pr_nd[c(1,1),], pr_nd3)

  #Check for factor level handling
  occ_covs$occ_fac <- factor(sample(c('a','b','c'),N,replace=T))

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)
  stateformulas <- c('~occ_fac','~1','~1')
  fm <- occuMulti(detformulas, stateformulas, data = umf)

  nd <- siteCovs(umf)[1:2,]
  pr_nd <- predict(fm, type='state', newdata=nd, se=F)$Predicted

  nd2 <- data.frame(occ_fac=factor(c('a','b'),levels=c('a','b','c')))
  pr_nd2 <- predict(fm, type='state', newdata=nd2, se=F)$Predicted

  expect_equivalent(pr_nd, pr_nd2[c(2,1),])

  nd3 <- data.frame(occ_fac=c('a','b'))
  pr_nd3 <- predict(fm, type='state', newdata=nd3, se=F)$Predicted

  expect_equivalent(pr_nd, pr_nd3[c(2,1),])

  nd4 <- data.frame(occ_fac=factor(c('a','d'),levels=c('a','d')))
  expect_error(predict(fm, type='state', newdata=nd4, se=F))

  #Check that predicting detection also works
  nd5 <- data.frame(det_cov1 = rnorm(5))
  pr_nd5 <- predict(fm, type='det', newdata=nd5)
  expect_equivalent(sapply(pr_nd5, nrow), c(5,5))
  expect_equivalent(pr_nd5$sp1$Predicted[1], 0.1680881)

  })
})

test_that("occuMulti penalized likelihood works",{

  set.seed(123)
  N <- 100; nspecies <- 3; J <- 5
  occ_covs <- as.data.frame(matrix(rnorm(N * 10),ncol=10))
  names(occ_covs) <- paste('occ_cov',1:10,sep='')

  det_covs <- list()
  for (i in 1:nspecies){
    det_covs[[i]] <- matrix(rnorm(N*J),nrow=N)
  }
  names(det_covs) <- paste('det_cov',1:nspecies,sep='')

  #True vals
  beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
  f1 <- beta[1] + beta[2]*occ_covs$occ_cov1
  f2 <- beta[3] + beta[4]*occ_covs$occ_cov2
  f3 <- beta[5] + beta[6]*occ_covs$occ_cov3
  f4 <- beta[7]
  f5 <- beta[8]
  f6 <- beta[9]
  f7 <- beta[10]
  f <- cbind(f1,f2,f3,f4,f5,f6,f7)
  z <- expand.grid(rep(list(1:0),nspecies))[,nspecies:1]
  colnames(z) <- paste('sp',1:nspecies,sep='')
  dm <- model.matrix(as.formula(paste0("~.^",nspecies,"-1")),z)
  psi <- exp(f %*% t(dm))
  psi <- psi/rowSums(psi)

  #True state
  ztruth <- matrix(NA,nrow=N,ncol=nspecies)
  for (i in 1:N){
    ztruth[i,] <- as.matrix(z[sample(8,1,prob=psi[i,]),])
  }
  p_true <- c(0.6,0.7,0.5)

  y <- list()
  for (i in 1:nspecies){
    y[[i]] <- matrix(NA,N,J)
    for (j in 1:N){
      for (k in 1:J){
        y[[i]][j,k] <- rbinom(1,1,ztruth[j,i]*p_true[i])
      }
    }
  }
  names(y) <- c('coyote','tiger','bear')

  umf = unmarkedFrameOccuMulti(y=y,siteCovs=occ_covs,obsCovs=det_covs)
  occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')
  detFormulas <- c('~1','~1','~1')

  nul <- capture_output({

  fm <- occuMulti(detFormulas,occFormulas,umf)
  fm_pen <- occuMulti(detFormulas, occFormulas, data = umf, penalty=1, boot=5)

  expect_equivalent(coef(fm_pen)[c(1,5)], c(0.5014605, -0.1078711), tol=1e-5)
  expect_equivalent(length(fm_pen@bootstrapSamples), 5)
  expect_equivalent(vcov(fm_pen), fm_pen@covMatBS)
  expect_equivalent(fm_pen@estimates@estimates$state@covMat,
                     fm_pen@estimates@estimates$state@covMatBS)

  set.seed(123)
  opt_fit  <- optimizePenalty(fm, penalties=c(0,1), boot=2)
  expect_equal(opt_fit@call$penalty, 1)

  })


})

test_that("Mismatched NAs are identified in unmarkedFrameOccuMulti",{
  y1 <- matrix(rbinom(10,1,0.5), nrow=5)
  y1[1,1] <- NA

  y2 <- matrix(rbinom(10,1,0.5), nrow=5)
  y2[1,1] <- NA

  y3 <- matrix(rbinom(10,1,0.5), nrow=5)
  y3[1,1] <- NA
  y3[5,1] <- NA

  ylist <- list(y1=y1,y2=y2,y3=y3)

  expect_warning(umf <- unmarkedFrameOccuMulti(y=ylist))

  pre_na <- sapply(ylist, function(x) sum(is.na(x)))
  post_na <- sapply(umf@ylist, function(x) sum(is.na(x)))

  expect_true(any(pre_na[1] != pre_na))
  expect_true(!any(post_na[1] != post_na))
})

test_that("R and C++ engines give same results",{

  y <- list(matrix(rep(0:1,10)[1:10],5,2),
            matrix(rep(0:1,10)[1:10],5,2))

  set.seed(123)
  N <- dim(y[[1]])[1]
  J <- dim(y[[1]])[2]
  occ_covs <- as.data.frame(matrix(rnorm(N * 3),ncol=3))
  names(occ_covs) <- paste('occ_cov',1:3,sep='')

  det_covs <- as.data.frame(matrix(rnorm(N*J*2),ncol=2))
  names(det_covs) <- paste('det_cov',1:2,sep='')

  umf <- unmarkedFrameOccuMulti(y = y, siteCovs = occ_covs, obsCovs = det_covs)
  stateformulas <- c('~occ_cov1','~occ_cov2','~occ_cov3')
  detformulas <- c('~det_cov1','~det_cov2')

  fm <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE,
                  control=list(maxit=1))
  fmR <- occuMulti(detformulas, stateformulas, data = umf, se=FALSE,
                  engine="R",control=list(maxit=1))
  expect_equal(coef(fm), coef(fmR))

})
