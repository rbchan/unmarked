
test.colext.na <- function() 
{

    y <- matrix(c(
        0,0, 1,1, 0,1,
        1,1, 0,0, 0,0,
        0,0, 0,0, 0,0), nrow=3, ncol=6, byrow=TRUE)
    
    y1 <- y
    y1[1,3] <- NA
    y1[2,5:6] <- NA
    
    umf1 <- unmarkedMultFrame(y=y1, numPrimary=3)
    fm1 <- colext(~1, ~1, ~1, ~1, umf1)
    
    
    oc <- y1 + -1:1
    umf2 <- unmarkedMultFrame(y=y, obsCovs=list(oc=oc), numPrimary=3)
    fm2 <- colext(~1, ~1, ~1, ~oc, umf2)


    ysc <- matrix(1:3, 3, 3, byrow=TRUE)
    ysc[1,1] <- NA    
    umf3 <- unmarkedMultFrame(y=y, yearlySiteCovs=list(ysc=ysc), numPrimary=3)
    checkException(fm3.1 <- colext(~1, ~1, ~ysc, ~1, umf3))
    checkException(fm3.2 <- colext(~1, ~ysc, ~1, ~1, umf3))

    ysc <- matrix(1:3, 3, 3, byrow=TRUE)
    ysc[1,1] <- NA
    y4 <- y
    y4[1,1:2] <- NA
    ysc4 <- ysc
    ysc4[1,1] <- 1    
    umf4 <- unmarkedMultFrame(y=y4, yearlySiteCovs=list(ysc=ysc4), numPrimary=3)
    fm4.1 <- colext(~1, ~1, ~ysc, ~1, umf4)
    fm4.2 <- colext(~1, ~ysc, ~1, ~1, umf4)


}