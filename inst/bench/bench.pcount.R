
library(unmarked)
library(rbenchmark)

set.seed(324534)
R <- 50
J <- 5
y <- matrix(NA, R, J)
z <- rpois(R, 1)
y[] <- rbinom(R*J, z, 0.5)

umf <- unmarkedFramePCount(y=y)


pcount(~1~1, umf, K=50)
pcountRcpp(~1~1, umf, K=50)

pcount.test <- benchmark(pcount(~1~1,umf,K=50),
                         pcountRcpp(~1~1,umf,K=50),
                         columns=c("test", "replications", "elapsed",
                                   "relative"),
                         replications=500)
pcount.test


