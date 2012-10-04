

test.umarkedMultFrame.subset <- function() {

    y <- matrix(1:27, 3)
    sc <- data.frame(x1 = 1:3)
    ysc <- list(x2 = matrix(1:9, 3))
    oc <- list(x3 = matrix(1:27, 3))

    umf1 <- unmarkedMultFrame(
        y = y,
        siteCovs = sc,
        yearlySiteCovs = ysc,
        obsCovs = oc,
        numPrimary = 3)

    dat <- as(umf1, "data.frame")

    umf1.obs1 <- umf1[,1]
    checkEquals(umf1.obs1@y, y[,1:3])
    checkEquals(umf1.obs1@siteCovs, sc)
    checkEqualsNumeric(unlist(umf1.obs1@obsCovs),
                       as.numeric(t(oc[[1]][,1:3])))
    checkEqualsNumeric(unlist(umf1.obs1@yearlySiteCovs), ysc[[1]][,1])
    checkEquals(umf1.obs1@numPrimary, 1)

    umf1.obs1and3 <- umf1[,c(1,3)]

    umf1.site1 <- umf1[1,]
    checkEquals(umf1.site1@y, y[1,, drop=FALSE])
    checkEquals(umf1.site1@siteCovs, sc[1,, drop=FALSE])
    checkEqualsNumeric(unlist(umf1.site1@obsCovs), oc$x3[1,])
    checkEqualsNumeric(unlist(umf1.site1@yearlySiteCovs),
        ysc$x2[1,, drop=FALSE])
    checkEquals(umf1.site1@numPrimary, 3)

    umf1.sites1and3 <- umf1[c(1,3),]

    }
