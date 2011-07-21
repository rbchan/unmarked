###################################################
### chunk number 1: 
###################################################
#line 6 "jss617.Rnw"
options(device.ask.default=FALSE)
options(width=70)
options(prompt="R> ")


###################################################
### chunk number 2: 
###################################################
#line 727 "jss617.Rnw"
library("unmarked")
data("mallard")


###################################################
### chunk number 3: 
###################################################
#line 737 "jss617.Rnw"
mallard.y[1:5,]


###################################################
### chunk number 4: 
###################################################
#line 744 "jss617.Rnw"
mallard.site[1:5,]


###################################################
### chunk number 5: 
###################################################
#line 757 "jss617.Rnw"
mallard.obs$ivel[1:5,]
mallard.obs$date[1:5,]


###################################################
### chunk number 6: 
###################################################
#line 770 "jss617.Rnw"
mallardUMF <- unmarkedFramePCount(y = mallard.y,
    siteCovs = mallard.site, obsCovs = mallard.obs)


###################################################
### chunk number 7: 
###################################################
#line 775 "jss617.Rnw"
summary(mallardUMF)


###################################################
### chunk number 8: 
###################################################
#line 799 "jss617.Rnw"
data("ovendata")


###################################################
### chunk number 9: 
###################################################
#line 809 "jss617.Rnw"
options(width=50)


###################################################
### chunk number 10: 
###################################################
#line 813 "jss617.Rnw"
ovenFrame <- unmarkedFrameMPois(ovendata.list$data, siteCovs = ovendata.list$covariates, type="removal")
ovenFrame[1:5,]


###################################################
### chunk number 11: 
###################################################
#line 818 "jss617.Rnw"
options(width=70)


###################################################
### chunk number 12: 
###################################################
#line 833 "jss617.Rnw"
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
#line 846 "jss617.Rnw"
y
umf <- unmarkedMultFrame(y=y, numPrimary=2)


###################################################
### chunk number 14: 
###################################################
#line 866 "jss617.Rnw"
ovenFrame[1:5, 1:3]


###################################################
### chunk number 15: 
###################################################
#line 876 "jss617.Rnw"
sc <- siteCovs(ovenFrame)
sc[,2:3] <- scale(sc[,2:3])
siteCovs(ovenFrame) <- sc


###################################################
### chunk number 16: 
###################################################
#line 911 "jss617.Rnw"
fm.mallard.1 <- pcount(~ date + I(date^2) ~ elev + forest, data = mallardUMF, K=50)
fm.mallard.1


###################################################
### chunk number 17: 
###################################################
#line 922 "jss617.Rnw"
fm.mallard.2 <- pcount(~ date ~ elev + forest, data = mallardUMF, K=50)
fm.mallard.2


###################################################
### chunk number 18: 
###################################################
#line 943 "jss617.Rnw"
fm.oven.1 <- multinomPois(~ 1 ~ ufp + trba, ovenFrame)
fm.oven.1


###################################################
### chunk number 19: 
###################################################
#line 962 "jss617.Rnw"
coef(fm.mallard.2)
coef(fm.mallard.2, type = "state")


###################################################
### chunk number 20: 
###################################################
#line 969 "jss617.Rnw"
names(fm.mallard.2)


###################################################
### chunk number 21: 
###################################################
#line 979 "jss617.Rnw"
vcov(fm.mallard.2, type = "det")


###################################################
### chunk number 22: 
###################################################
#line 986 "jss617.Rnw"
confint(fm.oven.1, type = "state", level = 0.95)


###################################################
### chunk number 23: 
###################################################
#line 995 "jss617.Rnw"
ci <- confint(fm.oven.1, type = "state", level = 0.95, method = "profile")


###################################################
### chunk number 24: 
###################################################
#line 998 "jss617.Rnw"
ci


###################################################
### chunk number 25: 
###################################################
#line 1019 "jss617.Rnw"
set.seed(1234)
fm.oven.1 <- nonparboot(fm.oven.1, B = 100)


###################################################
### chunk number 26: 
###################################################
#line 1024 "jss617.Rnw"
SE(fm.oven.1, type = "state")
SE(fm.oven.1, type = "state", method = "nonparboot")


###################################################
### chunk number 27: 
###################################################
#line 1032 "jss617.Rnw"
fm.oven.1 <- nonparboot(fm.oven.1, B = 100)


###################################################
### chunk number 28: 
###################################################
#line 1046 "jss617.Rnw"
(lc <- linearComb(fm.oven.1, type="state", coefficients=c(1,0.5,0)))


###################################################
### chunk number 29: 
###################################################
#line 1054 "jss617.Rnw"
(lc <- linearComb(fm.oven.1, type = "state",
    coefficients = matrix(c(1, 0.5, 0, 1, 1, 0), 2, 3, byrow = TRUE)))


###################################################
### chunk number 30: 
###################################################
#line 1063 "jss617.Rnw"
SE(lc, method = "nonparboot")


###################################################
### chunk number 31: 
###################################################
#line 1078 "jss617.Rnw"
(btlc <- backTransform(lc))
SE(btlc)
confint(btlc)


###################################################
### chunk number 32: 
###################################################
#line 1091 "jss617.Rnw"
fm.oven.2 <- update(fm.oven.1, formula = ~ 1 ~ ufp * trba)
fm.oven.3 <- update(fm.oven.1, formula = ~ 1 ~ ufp)
fm.oven.4 <- update(fm.oven.1, formula = ~ 1 ~ trba)
fm.oven.5 <- update(fm.oven.1, formula = ~ 1 ~ 1)


###################################################
### chunk number 33:  eval=FALSE
###################################################
## #line 1099 "jss617.Rnw"
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
#line 1119 "jss617.Rnw"
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
#line 1138 "jss617.Rnw"
options(width=60)


###################################################
### chunk number 36: 
###################################################
#line 1142 "jss617.Rnw"
fmList <- fitList("lam(ufp+trba)p(.)"=fm.oven.1, "lam(ufp*trba)p(.)"=fm.oven.2,
    "lam(ufp)p(.)"=fm.oven.3, "lam(trba)p(.)"=fm.oven.4, "lam(.)p(.)"=fm.oven.5)
modSel(fmList)


###################################################
### chunk number 37: 
###################################################
#line 1148 "jss617.Rnw"
options(width=70)


###################################################
### chunk number 38: 
###################################################
#line 1168 "jss617.Rnw"
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
#line 1181 "jss617.Rnw"
set.seed(1234)
plot(pb, main = "", xlab=expression(chi^2))


###################################################
### chunk number 40: 
###################################################
#line 1191 "jss617.Rnw"
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
## #line 1200 "jss617.Rnw"
## plot(pb, main = "")


