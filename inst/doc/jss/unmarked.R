###################################################
### chunk number 1: 
###################################################
#line 6 "unmarked.Rnw"
options(device.ask.default=FALSE)
options(width=70)
options(prompt="R> ")


###################################################
### chunk number 2: 
###################################################
#line 727 "unmarked.Rnw"
library("unmarked")
data("mallard")


###################################################
### chunk number 3: 
###################################################
#line 737 "unmarked.Rnw"
mallard.y[1:5,]


###################################################
### chunk number 4: 
###################################################
#line 744 "unmarked.Rnw"
mallard.site[1:5,]


###################################################
### chunk number 5: 
###################################################
#line 757 "unmarked.Rnw"
mallard.obs$ivel[1:5,]
mallard.obs$date[1:5,]


###################################################
### chunk number 6: 
###################################################
#line 770 "unmarked.Rnw"
mallardUMF <- unmarkedFramePCount(y = mallard.y,
    siteCovs = mallard.site, obsCovs = mallard.obs)


###################################################
### chunk number 7: 
###################################################
#line 775 "unmarked.Rnw"
summary(mallardUMF)


###################################################
### chunk number 8: 
###################################################
#line 799 "unmarked.Rnw"
data("ovendata")


###################################################
### chunk number 9: 
###################################################
#line 809 "unmarked.Rnw"
options(width=50)


###################################################
### chunk number 10: 
###################################################
#line 813 "unmarked.Rnw"
ovenFrame <- unmarkedFrameMPois(ovendata.list$data, siteCovs = ovendata.list$covariates, type="removal")
ovenFrame[1:5,]


###################################################
### chunk number 11: 
###################################################
#line 818 "unmarked.Rnw"
options(width=70)


###################################################
### chunk number 12: 
###################################################
#line 833 "unmarked.Rnw"
M <- 4
J <- 2
T <- 3
y <- matrix(0:1, nrow=M, ncol=J*T)
rownames(y) <- paste("site", 1:M, sep="")
cn <- paste("season", rep(1:2, each=3), sep="")
cn <- paste(cn, paste("v", rep(1:3, 2), sep=""), sep="")
colnames(y) <- cn


###################################################
### chunk number 13: 
###################################################
#line 846 "unmarked.Rnw"
y
umf <- unmarkedMultFrame(y=y, numPrimary=2)


###################################################
### chunk number 14: 
###################################################
#line 866 "unmarked.Rnw"
ovenFrame[1:5, 1:3]


###################################################
### chunk number 15: 
###################################################
#line 876 "unmarked.Rnw"
sc <- siteCovs(ovenFrame)
sc[,2:3] <- scale(sc[,2:3])
siteCovs(ovenFrame) <- sc


###################################################
### chunk number 16: 
###################################################
#line 911 "unmarked.Rnw"
fm.mallard.1 <- pcount(~ date + I(date^2) ~ elev + forest, data = mallardUMF, K=50)
fm.mallard.1


###################################################
### chunk number 17: 
###################################################
#line 922 "unmarked.Rnw"
fm.mallard.2 <- pcount(~ date ~ elev + forest, data = mallardUMF, K=50)
fm.mallard.2


###################################################
### chunk number 18: 
###################################################
#line 943 "unmarked.Rnw"
fm.oven.1 <- multinomPois(~ 1 ~ ufp + trba, ovenFrame)
fm.oven.1


###################################################
### chunk number 19: 
###################################################
#line 962 "unmarked.Rnw"
coef(fm.mallard.2)
coef(fm.mallard.2, type = "state")


###################################################
### chunk number 20: 
###################################################
#line 969 "unmarked.Rnw"
names(fm.mallard.2)


###################################################
### chunk number 21: 
###################################################
#line 979 "unmarked.Rnw"
vcov(fm.mallard.2, type = "det")


###################################################
### chunk number 22: 
###################################################
#line 986 "unmarked.Rnw"
confint(fm.oven.1, type = "state", level = 0.95)


###################################################
### chunk number 23: 
###################################################
#line 995 "unmarked.Rnw"
ci <- confint(fm.oven.1, type = "state", level = 0.95, method = "profile")


###################################################
### chunk number 24: 
###################################################
#line 998 "unmarked.Rnw"
ci


###################################################
### chunk number 25: 
###################################################
#line 1019 "unmarked.Rnw"
set.seed(1234)
fm.oven.1 <- nonparboot(fm.oven.1, B = 100)


###################################################
### chunk number 26: 
###################################################
#line 1024 "unmarked.Rnw"
SE(fm.oven.1, type = "state")
SE(fm.oven.1, type = "state", method = "nonparboot")


###################################################
### chunk number 27: 
###################################################
#line 1032 "unmarked.Rnw"
fm.oven.1 <- nonparboot(fm.oven.1, B = 100)


###################################################
### chunk number 28: 
###################################################
#line 1046 "unmarked.Rnw"
(lc <- linearComb(fm.oven.1, type="state", coefficients=c(1,0.5,0)))


###################################################
### chunk number 29: 
###################################################
#line 1054 "unmarked.Rnw"
(lc <- linearComb(fm.oven.1, type = "state",
    coefficients = matrix(c(1, 0.5, 0, 1, 1, 0), 2, 3, byrow = TRUE)))


###################################################
### chunk number 30: 
###################################################
#line 1063 "unmarked.Rnw"
SE(lc, method = "nonparboot")


###################################################
### chunk number 31: 
###################################################
#line 1078 "unmarked.Rnw"
(btlc <- backTransform(lc))
SE(btlc)
confint(btlc)


###################################################
### chunk number 32: 
###################################################
#line 1091 "unmarked.Rnw"
fm.oven.2 <- update(fm.oven.1, formula = ~ 1 ~ ufp * trba)
fm.oven.3 <- update(fm.oven.1, formula = ~ 1 ~ ufp)
fm.oven.4 <- update(fm.oven.1, formula = ~ 1 ~ trba)
fm.oven.5 <- update(fm.oven.1, formula = ~ 1 ~ 1)


###################################################
### chunk number 33:  eval=FALSE
###################################################
## #line 1099 "unmarked.Rnw"
## preddata <- predict(fm.oven.4, type = "state", appendData = TRUE)
## library(ggplot2)
## qplot(trba, Predicted, data = preddata, geom = "line",
##       xlab = "Scaled Basal Tree Area",
##       ylab = "Estimated Abundance") +
##   geom_ribbon(aes(x = trba, ymin = lower, ymax = upper),
##               alpha = 0.1) + theme_bw()


###################################################
### chunk number 34: 
###################################################
#line 1119 "unmarked.Rnw"
preddata <- predict(fm.oven.4, type = "state", appendData = TRUE)
library(ggplot2)
print(qplot(trba, Predicted, data = preddata, geom = "line",
      xlab = "Scaled Basal Tree Area",
      ylab = "Estimated Abundance") +
  geom_ribbon(aes(x = trba, ymin = lower, ymax = upper),
              alpha = 0.1) + theme_bw())


###################################################
### chunk number 35: 
###################################################
#line 1138 "unmarked.Rnw"
options(width=60)


###################################################
### chunk number 36: 
###################################################
#line 1142 "unmarked.Rnw"
fmList <- fitList("lam(ufp+trba)p(.)"=fm.oven.1, "lam(ufp*trba)p(.)"=fm.oven.2,
    "lam(ufp)p(.)"=fm.oven.3, "lam(trba)p(.)"=fm.oven.4, "lam(.)p(.)"=fm.oven.5)
modSel(fmList)


###################################################
### chunk number 37: 
###################################################
#line 1148 "unmarked.Rnw"
options(width=70)


###################################################
### chunk number 38: 
###################################################
#line 1168 "unmarked.Rnw"
set.seed(1234)
chisq <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((observed - expected)^2 / expected)
    }
pb <- parboot(fm.oven.1, statistic = chisq, nsim = 200)


###################################################
### chunk number 39: 
###################################################
#line 1181 "unmarked.Rnw"
set.seed(1234)
plot(pb, main = "", xlab=expression(chi^2))


###################################################
### chunk number 40: 
###################################################
#line 1191 "unmarked.Rnw"
set.seed(1234)
chisq <- function(fm) {
    observed <- getY(fm@data)
    expected <- fitted(fm)
    sum((observed - expected)^2 / expected)
    }
pb <- parboot(fm.oven.1, statistic = chisq, nsim = 200)


###################################################
### chunk number 41:  eval=FALSE
###################################################
## #line 1200 "unmarked.Rnw"
## plot(pb, main = "")


