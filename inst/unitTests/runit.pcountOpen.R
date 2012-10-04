

test.pcountOpen.null <- function()
{
  y <- matrix(c(
      3, 2, 1, 4,
      3, 4, 2, 1,
      0, 1, 2, 3,
      5, 3, 3, 4,
      2, 4, 3, 3), 5, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1))
  obsCovs <- data.frame(o1 = 1:20)
  umf <- unmarkedFramePCO(y = y, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, K=10,
      starts=c(1, 0, 0, 7))
  checkEqualsNumeric(coef(fm1), c(0.9565118, 0.2743964, 0.1349845, 7.0041091),
      tol = 1e-5)

  fm2 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="gamma", K=10)
  checkEqualsNumeric(coef(fm2), c(1.821936, 8.741870, -0.287360),
      tol = 1e-3)

  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf, se=FALSE, fix="omega", K=10)
  checkEqualsNumeric(coef(fm3), c(1.8111032, -0.6342198, -0.4520630),
      tol = 1e-5)

}






test.pcountOpen.na <- function()
{
  y1 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA,
      2, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf1 <- unmarkedFramePCO(y = y1, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm1 <- pcountOpen(~1, ~1, ~1, ~1, data = umf1, se=FALSE, K=10,
      starts=c(1.6, 0.24, 1.16, -0.268))
  checkEqualsNumeric(coef(fm1),
      c(1.49434036, 0.44381407, 0.80682012, 0.06490056), tol = 1e-5)

  y2 <- matrix(c(
      1, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1,
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  oc <- y1 + -2:3

  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- list(o1 = oc)
  ysc <- list(o2=oc)
  umf2 <- unmarkedFramePCO(y = y2, siteCovs = siteCovs, obsCovs = obsCovs,
      yearlySiteCovs=ysc, numPrimary=4)

  fm2.1 <- pcountOpen(~1, ~1, ~1, ~o1, data = umf2, se=FALSE, K=10,
      starts=c(1.4, -1.3, 1.8, -1.1, 0.7))
  checkEqualsNumeric(coef(fm2.1),
      c(1.2957439, -8.3373450, 2.2840248, -0.6967546, 1.1605447), tol = 1e-4)

  fm2.2 <- pcountOpen(~1, ~1, ~o2, ~1, data = umf2, se=FALSE, K=10,
      starts=c(1.4, -1.3, 1.8, -1.1, 0.7))
  checkEqualsNumeric(coef(fm2.2),
      c(1.36621986, 0.88669259, -2.46690971, -8.93330624, 0.02535309),
      tol = 1e-5)

  fm2.3 <- pcountOpen(~1, ~o2, ~1, ~1, data = umf2, se=FALSE, K=10,
      starts=c(1, 0, 0, -5, -1))
  checkEqualsNumeric(coef(fm2.3),
      c(0.7038926, 0.5306853, -0.2339208, -1.8501411, 4.5668710),
                     tol = 1e-2)

  y3 <- matrix(c(
      NA, 2, 1, 4,
      3, NA, 2, 1,
      0, 1, 2, NA,
      5, NA, 3, NA,
      NA, NA, 3, NA,
      NA, NA, NA, NA), 6, 4, byrow=TRUE)
  siteCovs <- data.frame(x = c(0,2,3,4,1,1))
  obsCovs <- data.frame(o1 = 1:24)
  umf3 <- unmarkedFramePCO(y = y3, siteCovs = siteCovs, obsCovs = obsCovs,
      numPrimary=4)

  fm3 <- pcountOpen(~1, ~1, ~1, ~1, data = umf3, se=FALSE, K=10,
      starts=c(1.5, 0, 1, 0))
  checkEqualsNumeric(coef(fm3), c(1.4751002, 0.4217504, 0.7234106, 0.1834803),
      tol = 1e-5)
  checkEquals(fm3@sitesRemoved, 6)



  y4 <- matrix(c(
      NA, 2, 1, 4,
      3, 1, 2, 1,
      0, 1, 2, 1,
      5, 1, 3, 1,
      1, 1, 3, 1,
      2, 1, 1, 1), 6, 4, byrow=TRUE)
  go4 <- matrix(c(
      NA, NA, NA, # remove y[1, 2:4]
      NA, 1, 2,   # remove y[2, 2]
      0, NA, 2,   # remove y[3, 3]. Creates an interior NA
      5, 1, NA,   # remove y[4, 4]. Creates an end NA
      1, NA, NA,  # remove y[5, 3:4]
      NA, NA, 1), 6, 3, byrow=TRUE)
  o2y <- matrix(c(
      1, 0, 0,
      0, 1, 0,
      0, 0, 1), 3, 3, byrow=TRUE)
  y4.na <- is.na(go4) %*% o2y
  y4.2 <- y4
  y4.2[,-1][y4.na>0] <- NA
  y4.2

  umf4 <- unmarkedFramePCO(y=y4, yearlySiteCovs=list(go4=cbind(go4, 1)),
      numPrimary=4)

  fm4.1 <- pcountOpen(~1, ~go4, ~1, ~1, umf4, se=FALSE,
      starts=c(.8, .5, -.3, -1.5, 6))
  checkEquals(fm4.1@sitesRemoved, 1)

  fm4.2 <- pcountOpen(~1, ~1, ~go4, ~1, umf4, se=FALSE,
      starts=c(.8, 0, 5, -5, 7))
  checkEquals(fm4.2@sitesRemoved, 1)



    # Now with secondary sampling periods

    y5 <- matrix(c(
        2,2,  2,2,  1,0,
        3,2,  1,1,  2,2,
        0,0,  1,1,  1,1,
        3,3,  2,0,  3,3), 4, 6, byrow=TRUE)

    umf5 <- unmarkedFramePCO(y=y5, numPrimary=3)
    fm5 <- pcountOpen(~1, ~1, ~1, ~1, umf5, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm5),
        c(0.7269958, -0.3484145, 0.1494188, 1.9391898), tol=1e-5)

    y6 <- y5
    y6[1,1] <- y6[2,3:4] <- y6[3,5:6] <- y6[4,6] <- NA
    umf6 <- unmarkedFramePCO(y=y6, numPrimary=3)
    fm6 <- pcountOpen(~1, ~1, ~1, ~1, umf6, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm6),
        c(0.7945817, -0.4340502, 0.5614526, 1.4161393), tol=1e-5)

    y7 <- y5
    oc7 <- y6 + -2:1
    umf7 <- unmarkedFramePCO(y=y7, obsCovs=list(oc=oc7), numPrimary=3)
    fm7 <- pcountOpen(~1, ~1, ~1, ~oc, umf7, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm7),
        c(1.1985964, -8.9001805, 12.2205930, -0.8876678, 0.9525985), tol=1e-4)

    y8 <- y5
    ysc8 <- matrix(1:3, 4, 3, byrow=TRUE)
    ysc8[1,1] <- NA
    umf8 <- unmarkedFramePCO(y=y8, yearlySiteCovs=list(ysc=ysc8), numPrimary=3)
    fm8 <- pcountOpen(~1, ~1, ~ysc, ~1, umf8, se=FALSE, K=10)
    checkEqualsNumeric(coef(fm8),
        c(0.7362607, -0.4708421, -2.3317736, 1.7017999, 1.8114414), tol=1e-4)

}







test.pcountOpen.delta <- function()
{
    M <- 5
    T <- 4
    y <- matrix(c(
        NA, 2, 1, 4,
        3, NA, 2, 1,
        0, 1, NA, 3,
        5, 3, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    if(!exists("formatDelta"))
        formatDelta <- unmarked:::formatDelta
    dates <- matrix(c(1,3,5,7), M, T, byrow=TRUE)
    delta <- formatDelta(dates, is.na(y))
    ans <- matrix(c(
        1, 2, 2, 2,
        1, 2, 4, 2,
        1, 2, 2, 4,
        1, 2, 2, 2,
        1, 2, 2, 2), M, T, byrow=TRUE)

    checkEquals(delta, ans)

    dates2 <- matrix(c(
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8), M, T, byrow=TRUE)
    delta2 <- formatDelta(dates2, is.na(y))
    ans2 <- matrix(c(
        2, 3, 2, 2,
        1, 3, 5, 2,
        2, 2, 2, 4,
        1, 3, 2, 2,
        2, 3, 2, 2), M, T, byrow=TRUE)

    checkEquals(delta2, ans2)

    dates3 <- matrix(as.integer(c(
      2, NA, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8,
      1, 4, 6, 8,
      2, 4, 6, 8)), M, T, byrow=TRUE)
    checkException(unmarkedFramePCO(y=y, primaryPeriod=dates3, numPrimary=4))

    dates4 <- dates2
    dates4[is.na(y)] <- NA
    mode(dates4) <- "integer"
    delta4 <- formatDelta(dates4, is.na(y))
    umf <- unmarkedFramePCO(y=y, primaryPeriod=dates4, numPrimary=4)
    fm <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10, starts=c(1.2, 0, 1.4, 1.2))
    checkEqualsNumeric(coef(fm),
        c(1.2206233, -0.1280961, 0.5874789, 5.9916012), tol = 1e-5)

    y5 <- matrix(c(
        1, NA, 1, 4,
        NA, 3, 2, 1,
        0, 1, 2, NA,
        NA, NA, 3, NA,
        NA, 4, NA, NA), M, T, byrow=TRUE)
    dates5 <- matrix(c(
        2, NA, 6, 8,
        NA, 4, 6, 8,
        2, 4, 6, NA,
        NA, NA, 6, NA,
        2, 4, 6, 8), M, T, byrow=TRUE)
    ans5 <- matrix(c(
        1, NA, 4, 2,
        NA, 2, 2, 2,
        1, 2, 2, NA,
        NA, NA, 4, NA, # 4 not 5 b/c primary period 1 is day 2
        1, 2, 2, 2), M, T, byrow=TRUE)
    delta5 <- formatDelta(dates5, is.na(y5))
    checkEquals(delta5, ans5)

    dates6 <- y6 <- matrix(c(2L, 1L), 1, 2)
    checkException(unmarkedFramePCO(y=y6, primaryPeriod=dates6, numPrimary=2))



}






test.pcountOpen.secondSamps <- function()
{
    y <- matrix(c(
        0,0,  2,2,  3,2,  2,2,
        2,2,  2,1,  3,2,  1,1,
        1,0,  1,1,  0,0,  0,0,
        0,0,  0,0,  0,0,  0,0), nrow=4, ncol=8, byrow=TRUE)

    sc <- data.frame(x1 = 1:4, x2 = c('A','A','B','B'))

    oc <- list(
        x3 = matrix(1:8, nrow=4, ncol=8, byrow=TRUE),
        x4 = matrix(letters[1:8], nrow=4, ncol=8, byrow=TRUE))

    ysc <- list(
        x5 = matrix(c(
            1,2,3,4,
            1,2,3,4,
            1,2,3,4,
            1,2,3,4), nrow=4, ncol=4, byrow=TRUE))

    umf1 <- unmarkedFramePCO(y=y, siteCovs=sc, obsCovs=oc,
        yearlySiteCovs=ysc, numPrimary=4)

    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=10)
    checkEqualsNumeric(coef(m1),
        c(-0.2438797, -0.7838448, 0.5572557, 1.6925454), tol=1e-5)


    y2 <- y
    y2[1,1] <- NA
    umf2 <- unmarkedFramePCO(y=y2, siteCovs=sc, obsCovs=oc,
        yearlySiteCovs=ysc, numPrimary=4)

    m2 <- pcountOpen(~1, ~1, ~1, ~1, umf2, K=10)



}

