test.gdistsamp.covs <- function() {
    y <- matrix(rep(12:1, 10), 5, 6, byrow=TRUE)
    siteCovs <- data.frame(x = c(-2, 2, 3, 4, 1))
    umf <- unmarkedFrameGDS(y = y, siteCovs = siteCovs,
        dist.breaks=c(0, 5, 10, 15)/1000, survey="line", tlength=rep(1, 5),
        unitsIn="km", numPrimary=2)
    fm <- gdistsamp(~x, ~x, ~x, data = umf, starts=c(2,0,-1,0,2,-1),
                    K=50, se=FALSE)

    lam <- fm['lambda']
    phi <- fm['phi']
    det <- fm['det']

    checkEqualsNumeric(coef(lam), c(3.72714393, -0.04496395),
                       tolerance = 1e-4)
    checkEqualsNumeric(coef(phi), c(0.41446960, -0.06101506),
                       tolerance = 1e-4)
    checkEqualsNumeric(coef(det), c(10.138343, -3.724821),
                       tolerance = 1e-4)

    }


