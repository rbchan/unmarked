
# -------------------------- Null Poisson removal model ------------------------

set.seed(26)

n <- 50  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals availabe for detection
    
for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)
    Nleft1 <- N[i,] - y[i,,1]
    y[i,,2] <- rbinom(T, Nleft1, p)
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }
    
y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")


system.time(m1 <- gmultmix(~1, ~1, ~1, data=umf1)) #2.3

# Test 1
checkEqualsNumeric(coef(m1), c(1.3923561, -0.3183231, -0.7864098), 
    tolerance=1e-5)

SSE(m1)

(pb1 <- parboot(m1, nsim=50, report=5))
plot(pb1)




# -------------------------- Null NegBin removal model ------------------------

set.seed(73)

n <- 50  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3
alpha <- 2

y <- array(NA, c(n, T, J))
M <- rnbinom(n, mu=lam, size=alpha)   # Local population size
N <- matrix(NA, n, T)               # Individuals availabe for detection
    
for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)
    Nleft1 <- N[i,] - y[i,,1]
    y[i,,2] <- rbinom(T, Nleft1, p)
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }
    
y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])
umf2 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")

system.time(m2 <- gmultmix(~1, ~1, ~1, data=umf2, mixture="NB")) #2.3

backTransform(m2, type="alpha")

# Test
checkEqualsNumeric(coef(m2), c(1.118504, 1.414340, -1.394736, 1.056084),
    tol=1e-5)

(pb2 <- parboot(m2, nsim=50, report=5))
plot(pb2)








# --------------------- Poisson removal model w/ covariates -------------------

set.seed(37)

n <- 50   # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

sc <- rnorm(n)
ysc <- rnorm(n*T)
ysc <- matrix(ysc, n, T)
yr <- factor(rep(1:T, n))
oc <- rnorm(n*J*T)
oc <- array(oc, c(n, J, T))
int <- matrix(1:(T*J), nrow=n, ncol=T*J, byrow=TRUE)
pi <- array(NA, c(n, J, T))

lam <- exp(-1 + 1*sc)
phi <- plogis(2 + -2*ysc)
p <- plogis(1 + -1*oc)

y <- array(NA, c(n,J,T))
M <- rpois(n, lam)
N <- matrix(NA, n, T)

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi[i,])
    y[i,1,] <- rbinom(T, N[i,], p[i,1,])
    Nleft1 <- N[i,] - y[i,1,]
    y[i,2,] <- rbinom(T, Nleft1, p[i,2,])
    Nleft2 <- Nleft1 - y[i,2,]
    y[i,3,] <- rbinom(T, Nleft2, p[i,3,])
    }

umf3 <- unmarkedFrameGMM(y=matrix(y, nrow=n), 
    siteCovs = data.frame(sc=sc), 
    obsCovs=list(oc=matrix(oc, nrow=n), int=int),
    yearlySiteCovs=data.frame(ysc=as.numeric(t(ysc)), yr=yr), 
    numPrimary=T, type="removal")

(m3 <- gmultmix(~sc, ~ysc, ~oc, umf3))
#system.time(m3 <- gmultmix(~sc, ~ysc, ~oc, umf3)) # 4.8
            
# Test
checkEqualsNumeric(coef(m3), c(-1.2513974, 1.3585940, 2.2889517, -2.1197854, 
    1.0450782, -0.8627125), tol=1e-5)
    
(pb3 <- parboot(m3, nsim=50, report=5))    
    
     


umf4 <- unmarkedFrameGMM(y=matrix(y, nrow=n), 
    siteCovs = data.frame(sc=sc), 
    obsCovs=list(oc=matrix(oc, nrow=n), int=int),
    yearlySiteCovs=list(ysc=ysc), 
    numPrimary=T, type="removal")
