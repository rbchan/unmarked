
library(unmarked)
library(rbenchmark)

set.seed(324534)
R <- 500
J <- 10
y <- matrix(NA, R, J)
z <- rbinom(R, 1, 0.5)
y[] <- rbinom(R*J, z, 0.5)

umf <- unmarkedFrameOccu(y=y)


occu(~1~1, umf)
occuRcpp(~1~1, umf)

occu.test <- benchmark(occu(~1~1, umf), occuRcpp(~1~1, umf),
    columns=c("test", "replications", "elapsed", "relative"),
    replications=50)

occu.test


