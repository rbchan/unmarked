



test.gmultmix.fit.NA <- function() 
{
    y <- matrix(0:3, 5, 4)
    siteCovs <- data.frame(x = c(0,2,3,4,1))
    siteCovs[3,1] <- NA
    obsCovs <- data.frame(o1 = 1:20, o2 = exp(-5:4)/20)
    yrSiteCovs <- data.frame(yr=factor(rep(1:2, 5)))
    
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs, 
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm <- gmultmix(~x, ~yr, ~o1 + o2, data = umf)
    checkEquals(fm@sitesRemoved, 3)
    checkEqualsNumeric(coef(fm), 
        c(2.50638554, 0.06226627, 0.21787839, 6.46029769, -1.51885928, 
            -0.03409375, 0.43424295), tol = 1e-5)

    obsCovs[10,2] <- NA
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs, 
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm <- gmultmix(~x, ~yr, ~o1 + o2, data = umf)
    checkEquals(fm@sitesRemoved, 3)
    checkEqualsNumeric(coef(fm), 
        c(2.50638554, 0.06226627, 0.21787839, 6.46029769, -1.51885928, 
            -0.03409375, 0.43424295), tol = 1e-5)
            
    yrSiteCovs[2, 1] <- NA
    umf <- unmarkedFrameGMM(y = y, siteCovs = siteCovs, obsCovs = obsCovs, 
        yearlySiteCovs = yrSiteCovs, type="removal", numPrimary=2)
    fm <- gmultmix(~x, ~yr, ~o1 + o2, data = umf)
    checkEqualsNumeric(coef(fm), 
        c(1.17280104, 0.37694710, 2.38249795, 2.87354955, -0.83875134, 
            -0.08446507, 1.88056826), tol = 1e-5)


}
