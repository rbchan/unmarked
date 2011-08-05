
library(unmarked)
library(rbenchmark)

set.seed(324534)
R <- 500
J <- 10
y <- matrix(NA, R, J)
z <- rbinom(R, 1, 0.5)
y[] <- rbinom(R*J, z, 0.5)

umf <- unmarkedFrameOccu(y=y)


(fm1.1 <- occu(~1~1, umf))
(fm1.2 <- occuRcpp(~1~1, umf))

all.equal(fm1.1, fm1.2, tolerance=1e-6)
all.equal(coef(fm1.1), coef(fm1.2), tolerance=1e-6)
all.equal(vcov(fm1.1), vcov(fm1.2), tolerance=1e-6)

(fm1.1 <- occu(~1~1, umf, control=list(reltol=1e-50)))
(fm1.2 <- occuRcpp(~1~1, umf, control=list(reltol=1e-50)))

all.equal(fm1.1, fm1.2, tolerance=1e-5)
all.equal(coef(fm1.1), coef(fm1.2), tolerance=1e-6)


occu.test <- benchmark(occu(~1~1, umf), occuRcpp(~1~1, umf),
    columns=c("test", "replications", "elapsed", "relative"),
    replications=50)

occu.test




# With missing values

y2 <- y
y2[1,1:J] <- NA
y2[2,1] <- NA

umf2 <- unmarkedFrameOccu(y=y2)


(fm2.1 <- occu(~1~1, umf2))
(fm2.2 <- occuRcpp(~1~1, umf2))

all.equal(coef(fm2.1), coef(fm2.2), tolerance=1e-6)

occu.test2 <- benchmark(occu(~1~1, umf2), occuRcpp(~1~1, umf2),
    columns=c("test", "replications", "elapsed", "relative"),
    replications=50)
occu.test2






# With known occ


known <- rep(FALSE, nrow(y))
known[1:3] <- TRUE

(fm3.1 <- occu(~1~1, umf, knownOcc=known))
(fm3.2 <- occuRcpp(~1~1, umf, knownOcc=known))

all.equal(coef(fm3.1), coef(fm3.2), tolerance=1e-6)

occu.test2 <- benchmark(occu(~1~1, umf2), occuRcpp(~1~1, umf2),
    columns=c("test", "replications", "elapsed", "relative"),
    replications=50)
occu.test2

