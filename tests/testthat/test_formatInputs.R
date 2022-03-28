context("input formatting functions")

test_that("formatDistData function works",{
    dat <- data.frame(distance=1:100, site=gl(5, 20),
                      visit=factor(rep(1:4, each=5)))
    cutpt <- seq(0, 100, by=25)
    y <- formatDistData(dat, "distance", "site", cutpt)
    expect_equivalent(y, matrix(c(20,   0,   0,   0,
                                    5,  15,   0,   0,
                                    0,  10,  10,   0,
                                    0,   0,  15,   5,
                                    0,   0,   0,  20), 5, 4, byrow=TRUE))
    dat.bad <- dat
    dat.bad$distance <- as.character(dat$distance)
    expect_error(formatDistData(dat.bad, "distance", "site", cutpt))

    dat.bad <- dat
    dat.bad$site <- as.character(dat$site)
    y2 <- expect_warning(formatDistData(dat.bad, "distance", "site", cutpt))
    expect_equivalent(y2, matrix(c(20,   0,   0,   0,
                                    5,  15,   0,   0,
                                    0,  10,  10,   0,
                                    0,   0,  15,   5,
                                    0,   0,   0,  20), 5, 4, byrow=TRUE))

    y3 <- formatDistData(dat, "distance", "site", cutpt, "visit")
    expect_equivalent(y3, matrix(c(
5, 0, 0, 0,   5, 0, 0, 0,   5, 0, 0, 0,   5, 0, 0, 0,
5, 0, 0, 0,   0, 5, 0, 0,   0, 5, 0, 0,   0, 5, 0, 0,
0, 5, 0, 0,   0, 5, 0, 0,   0, 0, 5, 0,   0, 0, 5, 0,
0, 0, 5, 0,   0, 0, 5, 0,   0, 0, 5, 0,   0, 0, 0, 5,
0, 0, 0, 5,   0, 0, 0, 5,   0, 0, 0, 5,   0, 0, 0, 5), 5, 16, byrow=TRUE))

    effortMatrix <- matrix(ncol=4, nrow=5,c(1,0))
    y4 <- formatDistData(dat, "distance","site",cutpt, "visit",effortMatrix)
    expect_equivalent(y4, matrix(c(
      5, 0, 0, 0,   NA,NA,NA,NA,  5, 0, 0, 0,   NA,NA,NA,NA,
      NA,NA,NA,NA,  0, 5, 0, 0,   NA,NA,NA,NA,   0, 5, 0, 0,
      0, 5, 0, 0,   NA,NA,NA,NA,  0, 0, 5, 0,   NA,NA,NA,NA,
      NA,NA,NA,NA,   0, 0, 5, 0,   NA,NA,NA,NA,  0, 0, 0, 5,
      0, 0, 0, 5,   NA,NA,NA,NA,   0, 0, 0, 5,   NA,NA,NA,NA), 5, 16, byrow=TRUE))

    effortMatrix <- matrix(ncol=4, nrow=5,"a")
    expect_error(formatDistData(dat, "distance","site",cutpt, "visit",effortMatrix))
})

test_that("formatLong works correctly",{
  df <- read.csv(system.file("csv","frog2001pcru.csv", package = "unmarked"),
                 stringsAsFactors=TRUE)
  umf <- formatLong(df, type = "unmarkedFrameOccu")
  ## Add some assertions...

  # Try simple with dates
  test <- expand.grid(site = LETTERS[1:2], date = as.Date(c("2017-04-06", "2017-04-11")))
  test <- test[with(test, order(site, date)), ]
  set.seed(1231)
  test <- within(test, {
    # ocov = round(rnorm(nrow(test)), 2)
    y = rbinom(nrow(test), 1, 0.6)
  })
  withdate <- formatLong(test, type = "unmarkedFrameOccu")

  expect_equal(withdate,
              new("unmarkedFrameOccu", y = structure(c(1L, 0L, 1L, 1L), .Dim = c(2L, 2L)),
                  obsCovs = structure(list(JulianDate = structure(c(17262, 17267, 17262, 17267),
                                                                  class = "Date")),
                                      class = "data.frame",
                                      row.names = c(NA, -4L)),
                  siteCovs = NULL,
                  mapInfo = NULL,
                  obsToY = structure(c(1, 0, 0, 1), .Dim = c(2L, 2L))))

  test <- expand.grid(site = LETTERS[1:4], julian = c(13, 20, 26))
  test <- test[with(test, order(site, julian)), ]

  set.seed(42)
  test <- within(test, {
    obsfac = factor(sample(LETTERS[1:2], nrow(test), replace = TRUE))
    sitefac = factor(round(as.numeric(site)/5))
    ocov = round(rnorm(nrow(test)), 2)
    scov = 2 * as.numeric(test$site)
    y = rbinom(nrow(test), 1, 0.6)
  })

  withfac <- formatLong(test, type = "unmarkedFrameOccu")

  expect_equal(withfac,
              new("unmarkedFrameOccu",
                  y = structure(c(1L, 0L, 1L, 0L, 1L, 1L, 0L, 0L, 0L, 1L, 1L, 0L), .Dim = 4:3),
                  obsCovs = structure(list(ocov = c(1.51, -0.09, 2.02, -0.06, 1.3, 2.29, -1.39, -0.28,
                                                    -0.13, 0.64, -0.28, -2.66),
                                           obsfac = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L,
                                                                2L, 1L, 2L, 1L, 2L),
                                                              .Label = c("A", "B"), class = "factor"),
                                           JulianDate = c(13, 20, 26, 13, 20, 26, 13, 20, 26, 13, 20, 26)),
                                      class = "data.frame", row.names = c(NA, -12L)),
                  siteCovs = structure(list(scov = c(2, 4, 6, 8),
                                            sitefac = structure(c(1L, 1L, 2L, 2L),
                                                                .Label = c("0", "1"), class = "factor")),
                                       class = "data.frame", row.names = c(NA, -4L)),
                  mapInfo = NULL,
                  obsToY = structure(c(1, 0, 0, 0, 1, 0, 0, 0, 1), .Dim = c(3L, 3L))))

  # Compare manual and automatic unmarkedPCount frames
  # Manual creation from help
  R <- 4 # number of sites
  J <- 3 # number of visits
  y <- matrix(c(1,2,0,0,0,0,1,1,1,2,2,1), nrow=R, ncol=J, byrow=TRUE)
  site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
  obs.covs <- list(
    x3 = matrix(c(
      -1,0,1,
      -2,0,0,
      -3,1,0,
      0,0,0), nrow=R, ncol=J, byrow=TRUE),
    x4 = matrix(factor(c(
      'a','b','c',
      'd','b','a',
      'a','a','c',
      'a','b','a')), nrow=R, ncol=J, byrow=TRUE))
  expect_warning(umf <- unmarkedFramePCount(y=y, siteCovs=site.covs,
                             obsCovs=obs.covs))        # organize data
  # Corresponding long data.frame
  pcdf <- data.frame(site = rep(seq(R), each = J),
                     occasion = rep(1:J, R),
                     y = as.vector(t(y)),
                     x1 = rep(1:4, each = J),
                     x2 = factor(rep(c('A','B', 'A', 'B'), each = J)),
                     x3 = as.vector(t(obs.covs$x3)),
                     x4 = as.vector(t(obs.covs$x4)))
  expect_warning(umf1 <- formatLong(pcdf, type = "unmarkedFramePCount"))
  # formatLong tacks on JulianDate to obsCovs, so ignore this difference
  expect_equal(umf@y, umf1@y)
  expect_equal(umf@siteCovs, umf1@siteCovs)
  expect_equal(umf@obsCovs, umf1@obsCovs[, c("x3", "x4")])
  expect_equal(umf@obsToY, umf1@obsToY)

  # Compare manual and automatic open point count frame
  y1 <- matrix(c(
    0, 2, 3, 2, 0,
    2, 2, 3, 1, 1,
    1, 1, 0, 0, 3,
    0, 0, 0, 0, 0), nrow=4, ncol=5, byrow=TRUE)

  # Site-specific covariates
  sc1 <- data.frame(x1 = 1:4, x2 = factor(c('A','A','B','B')))

  # Observation-specific covariates
  oc1 <- list(
    x3 = matrix(1:5, nrow=4, ncol=5, byrow=TRUE),
    x4 = matrix(letters[1:5], nrow=4, ncol=5, byrow=TRUE))

  # Primary periods of surveys
  primaryPeriod1 <- matrix(as.integer(c(
    1, 2, 5, 7, 8,
    1, 2, 3, 4, 5,
    1, 2, 4, 5, 6,
    1, 3, 5, 6, 7)), nrow=4, ncol=5, byrow=TRUE)

  # Create the unmarkedFrame
  expect_warning(umf1 <- unmarkedFramePCO(y=y1, siteCovs=sc1, obsCovs=oc1, numPrimary=5,
                           primaryPeriod=primaryPeriod1))

  test <- data.frame(site = rep(1:4, each = 5),
                     obsnum = 1:5,
                     y = as.vector(t(y1)),
                     x1 = rep(1:4, each = 5),
                     x2 = factor(rep(c('A','A','B','B'), each = 5)),
                     x3 = 1:5,
                     x4 = letters[1:5])
  expect_warning(umf2 <- formatLong(test, type = "unmarkedFramePCO", numPrimary = 5,
                     primaryPeriod = primaryPeriod1))
  # formatLong tacks on JulianDate to obsCovs, so ignore this difference
  expect_equal(umf1@y, umf2@y)
  expect_equal(umf1@siteCovs, umf2@siteCovs)
  expect_equal(umf1@obsCovs, umf2@obsCovs[, c("x3", "x4")])
  expect_equal(umf1@obsToY, umf2@obsToY)
  expect_equal(umf1@primaryPeriod, umf2@primaryPeriod)

  # Compare manual and automatic unmarkedFrameDS object
  # Manual creation from help
  R <- 4 # number of sites
  J <- 3 # number of distance classes
  db <- c(0, 10, 20, 30) # distance break points
  y <- matrix(c(
    5,4,3, # 5 detections in 0-10 distance class at this transect
    0,0,0,
    2,1,1,
    1,1,0), nrow=R, ncol=J, byrow=TRUE)
  site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
  umf <- unmarkedFrameDS(y=y, siteCovs=site.covs, dist.breaks=db, survey="point",
                         unitsIn="m")            # organize data
  # Corresponding long data.frame
  dsdf <- data.frame(site = rep(seq(R), each = J),
                     occasion = rep(1:J, R),
                     y = as.vector(t(y)),
                     x1 = rep(1:4, each = J),
                     x2 = factor(rep(c('A','B', 'A', 'B'), each = J)))
  umf1 <- formatLong(dsdf, type = "unmarkedFrameDS", dist.breaks = db,
                     survey = "point", unitsIn = "m")
  expect_equal(umf, umf1)

})


test_that("formatMult works correctly",{
  test <- expand.grid(site = LETTERS[1:4], visit = 1:3, year = 2015:2016)
  test <- test[with(test, order(site, year, visit)), ]
  test <- test[, c("year", "site", "visit")]

  set.seed(18939)
  test <- within(test, {
    ocov = round(rnorm(nrow(test)), 2)
    scov = 2 * as.numeric(test$site)
    yscov = 1.3 * as.numeric(interaction(test$site, test$year))
    obsfac = factor(sample(LETTERS[1:2], nrow(test), replace = TRUE))
    sitefac = factor(round(as.numeric(site)/5))
    ysfac = factor(round(as.numeric(interaction(site, year))/10))
    y  = rpois(nrow(test), lambda = 2)
  })

  withfac <- formatMult(test)

  expect_equal(withfac,
              new("unmarkedMultFrame",
                  numPrimary = 2L,
                  yearlySiteCovs = structure(list(ysfac = structure(c(1L, 1L, 1L, 2L, 1L, 2L, 1L, 2L),
                                                                    .Label = c("0", "1"), class = "factor"),
                                                  yscov = c(1.3, 6.5, 2.6, 7.8, 3.9, 9.1, 5.2, 10.4)),
                                             class = "data.frame", row.names = c(NA, -8L)),
                  y = structure(c(0L, 1L, 3L, 3L, 2L, 2L, 2L, 1L, 1L, 0L, 3L, 3L, 1L, 1L, 0L, 0L, 3L, 1L, 2L,
                                  2L, 2L, 1L, 3L, 3L), .Dim = c(4L, 6L)),
                  obsCovs = structure(list(visit = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1,
                                                     2, 3, 1, 2, 3),
                                           obsfac = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L,
                                                                2L, 1L, 2L, 1L, 1L, 2L, 1L,
                                                                2L, 2L, 2L, 1L, 2L, 1L, 2L), .Label = c("A", "B"), class = "factor"),
                                           ocov = c(0.28, -1.41, -0.31, 0.05, -0.53, 0.84, -0.95, 1.63, 0.87,
                                                    1.03, 1.41, 1.25, -0.32, 0.11, -0.45, -0.83, 0.17, 0.28,
                                                    -0.13, -1.86, -1.82, 0.11, 1.29, -0.31)),
                                      class = "data.frame", row.names = c(NA, -24L)),
                  siteCovs = structure(list(sitefac = structure(c(1L, 1L, 2L, 2L),
                                                                .Label = c("0", "1"), class = "factor"),
                                            scov = c(2, 4, 6, 8)), class = "data.frame",
                                       row.names = c(NA, -4L)),
                  mapInfo = NULL,
                  obsToY = structure(c(1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                       0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                       0, 0, 0, 0, 1), .Dim = c(6L, 6L))))
})

test_that("csvToUMF function works",{


  expect_warning(umf <- csvToUMF(system.file("csv","csv_factor_test.csv",
                                  package = "unmarked"), type="unmarkedFrameOccu"))


  expect_warning(umf <- csvToUMF(system.file("csv","csv_factor_test.csv",
                  package = "unmarked"), type="unmarkedFrameOccu"))

  expect_equal(sapply(siteCovs(umf), class), c(elev="numeric", forest="factor"))
  expect_equal(sapply(obsCovs(umf), class), c(wind="numeric", rain="factor"))

  df <- as(umf, "data.frame")
  expect_equivalent(dim(df), c(20,11))
})
