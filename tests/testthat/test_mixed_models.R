context("mixed model tools")
skip_on_cran()

test_that("get_reTrms matches lme4::mkReTrms", {

  #skip_if(!requireNamespace("lme4", quietly=TRUE), 
  #        "lme4 package unavailable")
  set.seed(123) 
  dat <- data.frame(x = rnorm(20), y = rnorm(20), z = factor(sample(letters[1:3], 20, replace=T)),
                    group = factor(sample(letters[4:6], 20, replace=T)),
                    id = factor(sample(letters[7:9], 20, replace=T)))
  
  load('lme4_output.Rdata')
  form1 <- ~x + (1|group)
  #l1 <- lme4::mkReTrms(lme4::findbars(form1), dat)
  r1 <- get_reTrms(form1, dat)
  expect_identical(r1$Z, Matrix::t(l1$Zt))
  expect_identical(r1$cnms, l1$cnms)
  attributes(l1$flist)$assign <- NULL
  expect_identical(r1$flist, l1$flist)

  form2 <- ~x + (x||group)
  #l2 <- lme4::mkReTrms(lme4::findbars(form2), dat)
  r2 <- get_reTrms(form2, dat)
  expect_identical(r2$Z, Matrix::t(l2$Zt))
  expect_identical(r2$cnms, l2$cnms)
  attributes(l2$flist)$assign <- NULL
  expect_identical(r2$flist, l2$flist)

  form3 <- ~x + (x||group) + (1|id)
  #l3 <- lme4::mkReTrms(lme4::findbars(form3), dat)
  r3 <- get_reTrms(form3, dat)
  expect_identical(r3$Z, Matrix::t(l3$Zt))
  expect_identical(r3$cnms, l3$cnms)
  attributes(l3$flist)$assign <- NULL
  expect_identical(r3$flist, l3$flist)

  form4 <- ~x + (x||group) + (y||id)
  #l4 <- lme4::mkReTrms(lme4::findbars(form4), dat)
  r4 <- get_reTrms(form4, dat)
  expect_identical(r4$Z, Matrix::t(l4$Zt))
  expect_identical(r4$cnms, l4$cnms)
  attributes(l4$flist)$assign <- NULL
  expect_identical(r4$flist, l4$flist)

  form5 <- ~x + (x*y||group) + (y||id)
  #l5 <- lme4::mkReTrms(lme4::findbars(form5), dat)
  r5 <- get_reTrms(form5, dat)
  expect_identical(r5$Z, Matrix::t(l5$Zt))
  expect_identical(r5$cnms, l5$cnms)
  attributes(l5$flist)$assign <- NULL
  expect_identical(r5$flist, l5$flist)

  form6 <- ~(1|group)
  #l6 <- lme4::mkReTrms(lme4::findbars(form6), dat)
  r6 <- get_reTrms(form6, dat)
  expect_identical(r6$Z, Matrix::t(l6$Zt))
  expect_identical(r6$cnms, l6$cnms)
  attributes(l6$flist)$assign <- NULL
  expect_identical(r6$flist, l6$flist)

  form7 <- ~(1|group) + x
  #l7 <- lme4::mkReTrms(lme4::findbars(form7), dat)
  r7 <- get_reTrms(form7, dat)
  expect_identical(r7$Z, Matrix::t(l7$Zt))
  expect_identical(r7$cnms, l7$cnms)
  attributes(l7$flist)$assign <- NULL
  expect_identical(r7$flist, l7$flist)

  # Check that unused factor levels aren't dropped
  dat2 <- data.frame(x = c(0.1, 0.2, -0.1),
                    group = factor(c("a","b","c"), levels=c("a","b","c","d")))
  form8 <- ~x + (1|group)
  #l8 <- lme4::mkReTrms(lme4::findbars(form8), dat2, drop.unused.levels=FALSE)
  r8 <- get_reTrms(form8, dat2)
  expect_identical(r8$Z, Matrix::t(l8$Zt))

  # Check that get_Z handles newdata
  form9 <- ~x + (x||group) + (1|id)
  nd <- data.frame(x=c(0.5,1), group=c("e","d"), id=c("i","i"))
  Z <- get_Z(form9, dat, newdata=nd)
  expect_equivalent(as.matrix(Z), matrix(c(0,1,0,0,0.5,0,0,0,1,
                                      1,0,0,1,0,0,0,0,1), nrow=2, byrow=T))

  # New level
  nd <- data.frame(x=c(0.5,1), group=c("a","d"), id=c("i","i"))
  expect_error(get_Z(form9, dat, newdata=nd))
  
  #save(l1,l2,l3,l4,l5,l6,l7,l8, file='lme4_output.Rdata')
})

test_that("get_reTrms errors correctly", {
  set.seed(123) 
  dat <- data.frame(x = rnorm(20), y = rnorm(20), z = factor(sample(letters[1:3], 20, replace=T)),
                    group = factor(sample(letters[4:6], 20, replace=T)),
                    id = factor(sample(letters[7:9], 20, replace=T)))
  
  form1 <- ~x + (x+z||group) + (y||id)
  expect_error(get_reTrms(form1, dat))

  form1 <- ~x + (x|group:id)
  expect_error(get_reTrms(form1, dat))

  form1 <- ~x + (x|group/id)
  expect_error(get_reTrms(form1, dat))

  expect_identical(find_bars(NULL), NULL)
})
