test.emptyframe <- function() {
  checkException(umf <- unmarkedFrame())
}

test.onecol <- function() {
    checkException(umf <- unmarkedFrame(y=matrix(1)))
    }

test.frame <- function() {
  M <- 10
  J <- 3
  y <- matrix(rbinom(J * M, 1, 0.5), M, J)
  siteCovs <- data.frame(a = rnorm(M), b = factor(gl(2,5)))
  umf <- unmarkedFrame(y = y, siteCovs = siteCovs)
}


test.umfDS.args <- function() {
    y <- matrix(1, 1, 2)
    s <- "point"
    d <- 0:2
    uin <- "m"
    sc <- data.frame(1)
    oc <- matrix(1, 2)
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, obsCovs=oc, 
        survey=s, dist.breaks=d, unitsIn=uin))
    umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, dist.breaks=d, 
        unitsIn=uin)    
    checkException(obsCovs(umf) <- oc)
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, 
        dist.breaks=d))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, dist.breaks=d, 
        unitsIn=uin))    
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, 
        unitsIn=uin))    
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, 
        dist.breaks=0:3, unitsIn=uin))
    checkException(umf <- unmarkedFrameDS(y=y, siteCovs=sc, survey=s, 
        dist.breaks=1:3, unitsIn=uin))
             
    }    