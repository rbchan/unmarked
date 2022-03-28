context("simulate method")

test_that("simulate can generate new datasets from scratch",{

  set.seed(123)
  forms <- list(state=~elev, det=~1)
  design <- list(M=300, J=5)

  # Should write a better handler for this situation
  bad_forms <- list(occu=~elev, det=~1)
  expect_error(simulate("occu", formulas=bad_forms, design=design))

  # When no coefficients list provided
  nul <- capture_output(expect_error(simulate("occu", formulas=forms, design=design)))

  cf <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))
  umf <- simulate("occu", formulas=forms, design=design, coefs=cf)
  expect_is(umf, "unmarkedFrame")
  expect_equivalent(dim(umf@y), c(300,5))
  expect_equal(names(umf@siteCovs), "elev")

  fm <- occu(~1~elev, umf)
  expect_equivalent(coef(fm), c(-0.06492,-0.43037,0.0527354), tol=1e-4)

  # With guide
  set.seed(123)
  guide <- list(elev=list(dist=rnorm, mean=2, sd=0.5),
                landcover=factor(levels=c("forest","grass")))
  forms$state <- ~elev+landcover
  cf$state <- c(intercept=0, elev=-0.4, landcovergrass=0.5)
  umf2 <- simulate("occu", formulas=forms, design=design, coefs=cf, guide=guide)
  expect_equal(names(umf2@siteCovs), c("elev","landcover"))
  expect_true(is.factor(umf2@siteCovs$landcover))
  expect_equivalent(mean(umf2@siteCovs$elev), 2.01722, tol=1e-5)

  # pcount
  set.seed(123)
  cf$alpha <- c(alpha=0.5)
  umf3 <- simulate("pcount", formulas=forms, design=design, coefs=cf, guide=guide,
                   mixture="NB", K=10)
  fm2 <- pcount(~1~elev, umf3, mixture="NB", K=10)
  expect_equivalent(coef(fm2), c(-0.1775, -0.2528, -0.083, 0.5293), tol=1e-3)

  # distsamp
  set.seed(123)
  cf$alpha <- NULL
  cf$det[1] <- log(30)
  cf$state <- c(intercept=2, elev=0.5)
  forms$state <- ~elev
  umf4 <- simulate("distsamp", formulas=forms, design=design, coefs=cf,
                    dist.breaks=c(0,10,20,30,40,50), survey='point', unitsIn='m')
  fm <- distsamp(~1~elev, umf4)
  expect_equivalent(coef(fm), c(1.9389,0.5344,3.4521), tol=1e-4)

  # Mpois
  set.seed(123)
  cf$dist[1] <- 0
  cf$state <- c(intercept=1, elev=0.5)
  umf5 <- simulate("multinomPois", formulas=forms, design=design, coefs=cf,
                   guide=guide)
  fm <- multinomPois(~1~elev, umf5)
  expect_equivalent(coef(fm), c(0.98163,0.50477,3.3633), tol=1e-3)

  #colext
  set.seed(123)
  forms_colext <- list(psi=~elev, col=~1, ext=~1, det=~1)
  cf_colext <- list(psi=c(intercept=0, elev=0.5), col=c(intercept=0),
                    ext=c(intercept=0), det=c(intercept=0))
  design_colext <- list(M=300,T=3,J=5)
  umf6 <- simulate("colext", formulas=forms_colext, design=design_colext,
                   coefs=cf_colext)
  fm <- colext(~elev, ~1, ~1, ~1, umf6)
  expect_equivalent(coef(fm), c(0.1598,0.6468,-0.0097,-0.01665,-0.0104),
                     tol=1e-3)

  #occuTTD
  set.seed(123)
  cf_ttd <- cf_colext
  cf_ttd$det <- c(intercept=log(0.5))
  umf7 <- simulate("occuTTD", formulas=forms_colext, design=design_colext,
                   coefs=cf_ttd, surveyLength=3)
  fm <- occuTTD(~elev, ~1, ~1, ~1, umf7)
  expect_equivalent(coef(fm), c(-0.0434,0.5743,-0.0187,0.115,-0.672),
                     tol=1e-3)

  #gdistsamp
  set.seed(123)
  cf_gds <- list(det=c(intercept=log(30)), lambda=c(intercept=2, elev=0.5),
                 phi=c(intercept=0))
  forms_gds <- list(lambda=~elev, phi=~1, det=~1)
  umf8 <- simulate("gdistsamp", formulas=forms_gds, design=design_colext, coefs=cf_gds,
                    dist.breaks=c(0,10,20,30,40,50), survey='line',
                    tlength=rep(100,300), unitsIn='m')
  fm <- gdistsamp(~elev,~1,~1, umf8)
  expect_equivalent(coef(fm), c(1.98053,0.5268,-0.05892,3.4113), tol=1e-3)

  #gmultmix
  set.seed(123)
  cf_gmm <- list(det=c(intercept=0), lambda=c(intercept=2, elev=0.5),
                 phi=c(intercept=0))
  forms_gmm <- list(lambda=~elev, phi=~1, det=~1)
  umf9 <- simulate("gmultmix", formulas=forms_gmm, design=design_colext, coefs=cf_gmm,
                   type='removal')
  fm <- gmultmix(~elev,~1,~1, umf9)
  expect_equivalent(coef(fm), c(1.9529,0.5321,0.0529,-0.0373), tol=1e-4)

  #gpcount
  set.seed(123)
  umf10 <- simulate("gpcount", formulas=forms_gmm, design=list(M=50,J=5,T=3), coefs=cf_gmm,
                    K=10)
  fm <- gpcount(~elev,~1,~1, umf10, K=10)
  expect_equivalent(coef(fm), c(1.4994,0.4024,1.1351,0.0978), tol=1e-4)

  #pcountOpen
  set.seed(123)
  cf_pco <- list(lambda=c(intercept=2, elev=0.5), det=c(intercept=0),
                 gamma=c(intercept=0), omega=c(intercept=0))
  design_pco <- list(M=100,J=5,T=3)
  forms_pco <- list(lambda=~elev, det=~1, gamma=~1, omega=~1)
  umf11 <- simulate("pcountOpen", formulas=forms_pco, design=design_pco,
                    coefs=cf_pco, K=15)
  fm <- pcountOpen(~elev, ~1, ~1, ~1, data=umf11, K=15)
  expect_equivalent(coef(fm), c(1.7703,0.0427,-0.2768,0.1288,0.0245), tol=1e-4)

  #multmixOpen
  set.seed(123)
  umf12 <- simulate("multmixOpen", formulas=forms_pco, design=design_pco,
                    coefs=cf_pco, K=15, type='removal')
  expect_is(umf12, "unmarkedFrameMMO")
  #fm <- multmixOpen(~elev,~1,~1,~1, data=umf12, K=15)
  #expect_equivalent(coef(fm), c(1.8128,0.0171,-0.4220,0.1921,-0.1122),tol=1e-4)

  #distsampOpen
  set.seed(123)
  cf_dso <- cf_pco
  cf_pco$det <- c(intercept=log(30))
  design_dso <- list(M=200, J=5, T=5)
  umf13 <- simulate("distsampOpen", formulas=forms_pco, design=design_dso,
                    coefs=cf_dso, K=20, unitsIn='m',
                    survey='point', dist.breaks=c(0,10,20,30,40,50))
  expect_is(umf13, "unmarkedFrameDSO")
  #fm <- distsampOpen(~elev,~1,~1,~1, data=umf13, K=20)
  #expect_equivalent(coef(fm), c(1.70195,0.00067,-0.1141,0.09816,3.4179), tol=1e-4)

  # occuMulti
  set.seed(123)
  occFormulas <- c('~occ_cov1','~occ_cov2','~occ_cov3','~1','~1','~1','~1')
  detFormulas <- c('~1','~1','~1')
  beta <- c(0.5,0.2,0.4,0.5,-0.1,-0.3,0.2,0.1,-1,0.1)
  p_true <- c(0.6,0.7,0.5)

  cf <- list(state=beta, det=log(p_true/(1-p_true)))
  names(cf$state) <-  c("[sp1] intercept", "[sp1] occ_cov1",
                      "[sp2] intercept", "[sp2] occ_cov2",
                      "[sp3] intercept", "[sp3] occ_cov3",
                      "[sp1:sp2] intercept","[sp1:sp3] intercept",
                      "[sp2:sp3] intercept","[sp1:sp2:sp3] intercept")
  names(cf$det) <- c("[sp1] intercept", "[sp2] intercept", "[sp3] intercept")

  umf14 <- simulate("occuMulti", formulas=list(state=occFormulas, det=detFormulas),
                design=list(M=200, J=5), coefs=cf)
  fm <- occuMulti(detFormulas, occFormulas, umf14)
  expect_equivalent(coef(fm, 'det'), c(0.3650,0.8762,-0.04653), tol=1e-4)

  # occuMS
  set.seed(123)
  bstate <- c(-0.5, 1, -0.6, -0.7)
  bdet <- c(-0.4, 0, -1.09, -0.84)
  detformulas <- c('~V1','~1','~1')
  stateformulas <- c('~V1','~V2')
  forms <- list(det=detformulas, state=stateformulas)
  cf <- list(state=bstate, det=bdet)
  expect_warning(umf15 <- simulate("occuMS", formulas=forms, coefs=cf, design=list(M=500, J=5, T=1)))
  fm <- occuMS(forms$det, forms$state, data=umf15, parameterization="multinomial")
  expect_equivalent(coef(fm, 'state'), c(-0.437,0.767,-0.671,-0.595), tol=1e-3)

  # gdistremoval
  set.seed(123)
  formulas <- list(lambda=~sc1, rem=~oc1, dist=~1, phi=~1)
  cf <- list(lambda=c(intercept=log(5), sc1=0.7), dist=c(intercept=log(50)),
           rem=c(intercept=log(0.2/(1-0.2)), oc1=0.4))
  design <- list(M=500, Jdist=4, Jrem=5, T=1)
  umf16 <- simulate("gdistremoval", design=design, formulas=formulas, coefs=cf,
                 dist.breaks=c(0,25,50,75,100), unitsIn='m', output='abund',K=15)
  fm <- gdistremoval(~sc1, removalformula=~oc1, distanceformula=~1,
                     data=umf16,K=15)
  expect_is(fm, "unmarkedFitGDS")

})
