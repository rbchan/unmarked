
library(unmarked)

## Simulate no covariates, constant sampling period intervals,
## no secondary samples

sim1 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}





set.seed(3223)
nsim1 <- 50
simout1 <- matrix(NA, nsim1, 4)
colnames(simout1) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim1) {
    cat("sim", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    y.sim1 <- sim1(lambda, gamma, omega, p)
    umf1 <- unmarkedFramePCO(y = y.sim1, numPrimary=5)
    m1 <- pcountOpen(~1, ~1, ~1, ~1, umf1, K=15,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)),
        se=FALSE)
    e <- coef(m1)
    simout1[i, 1:2] <- exp(e[1:2])
    simout1[i, 3:4] <- plogis(e[3:4])
    cat("  beta.hat =", simout1[i,], "\n")
    }

#png("pcountOpenSim1.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout1[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout1[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout1[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout1[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
#dev.off()





## Simulate covariate model with constant intervals



sim2 <- function(lam=c(0,1), gam=c(-1,-1), om=c(2,-1), p=c(-1,1), M=100,
                 T=5)
{
    y <- gamma <- omega <- det <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    veght <- rnorm(M)
    isolation <- matrix(rnorm(M*T), M, T)
    time <- matrix(rnorm(M*T, 1), M, T)
    lambda <- exp(lam[1] + lam[2]*veght)
    gamma[] <- exp(gam[1] + gam[2]*isolation)
    omega[] <- plogis(om[1] + om[2]*isolation)
    det[] <- plogis(p[1] + p[2]*time)

    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega[,t])
        G[,t] <- rpois(M, gamma[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, det)
    return(list(y=y, covs=data.frame(veght=veght,
        isolation=isolation, time=time)))
}





nsim2 <- 10
simout2 <- matrix(NA, nsim2, 8)
colnames(simout2) <- c('lam0', 'lam1', 'gam0', 'gam1', 'om0', 'om1',
                       'p0', 'p1')
for(i in 1:nsim2) {
    cat("sim", i, "\n"); flush.console()
    lam <- c(-2, 1)
    gam <- c(-1, -1)
    om <- c(0, -1)
    p <- c(-1, 1)
    T <- 5
    sim2out <- sim2(lam, gam, om, p, T=T)
    y.sim2 <- sim2out$y
    covs <- sim2out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    yearlySiteCovs <- list(isolation = covs[,grep("isolation", cn)])
    obsCovs <- list(time = covs[,grep("time", cn)])
    umf2 <- unmarkedFramePCO(y = y.sim2, siteCovs=siteCovs,
        yearlySiteCovs=yearlySiteCovs, obsCovs=obsCovs, numPrimary=T)
    m2 <- pcountOpen(~veght, ~isolation, ~isolation, ~time, umf2,
                     K=30, se=F, starts=c(lam, gam, om, p))
    e <- coef(m2)
    simout2[i, ] <- e
    cat("  beta.hat =", e, "\n")
    }

#png("pcountOpenSim2.png", width=6, height=8, units="in", res=360)
par(mfrow=c(4,2))
hist(simout2[,1], xlab=expression(lambda)); abline(v=lam[1], lwd=2, col=4)
hist(simout2[,2], xlab=expression(lambda)); abline(v=lam[2], lwd=2, col=4)
hist(simout2[,3], xlab=expression(gamma)); abline(v=gam[1], lwd=2, col=4)
hist(simout2[,4], xlab=expression(gamma)); abline(v=gam[2], lwd=2, col=4)
hist(simout2[,5], xlab=expression(omega)); abline(v=om[1], lwd=2, col=4)
hist(simout2[,6], xlab=expression(omega)); abline(v=om[2], lwd=2, col=4)
hist(simout2[,7], xlab=expression(p)); abline(v=p[1], lwd=2, col=4)
hist(simout2[,8], xlab=expression(p)); abline(v=p[2], lwd=2, col=4)
#dev.off()













## Simulate uneven sampling period intervals with all dates[i,1]==1
set.seed(333)

sim3 <- function(lambda=4, gamma=0.1, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- date <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    date[,1] <- 1

    for(i in 1:M) {
    for(t in 2:T) {
        delta <- max(rpois(1, 3), 1)
        date[i, t] <- date[i, t-1] + delta
        S[i, t-1] <- rbinom(1, N[i, t-1], omega)
        G[i, t-1] <- rpois(1, gamma)
        N[i, t] <- S[i, t-1] + G[i, t-1]
        if(delta > 1) {
            for(d in 2:delta) {
                S[i, t-1] <- rbinom(1, N[i, t], omega)
                G[i, t-1] <- rpois(1, gamma)
                N[i, t] <- S[i, t-1] + G[i, t-1]
                }
            }
        }}
    y[] <- rbinom(M*T, N, p)
    mode(date) <- "integer"
    return(list(y=y, dates=date))
}










set.seed(373)
nsim3 <- 500
simout3 <- matrix(NA, nsim3, 4)
colnames(simout3) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim3) {
    cat("sim3", i, "\n"); flush.console()
    lambda <- 4
    gamma <- 0.3
    omega <- 0.7
    p <- 0.7
    T <- 5
    yd <- sim3(lambda, gamma, omega, p, M=100, T=5)
    y.sim3 <- yd$y
    dates3 <- yd$dates
    umf3 <- unmarkedFramePCO(y = y.sim3, primaryPeriod=dates3, numPrimary=T)
    m3 <- pcountOpen(~1, ~1, ~1, ~1, umf3, K=20,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m3)
    simout3[i, 1:2] <- exp(e[1:2])
    simout3[i, 3:4] <- plogis(e[3:4])
    }



png("pcountOpenSim3.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout3[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout3[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout3[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout3[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()














# Auto-regressive model




sim4 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma*N[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}



set.seed(3223)
nsim4 <- 500
simout4 <- matrix(NA, nsim4, 4)
colnames(simout4) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim4) {
    cat("sim4", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.7
    p <- 0.7
    T <- 5
    y.sim4 <- sim4(lambda, gamma, omega, p, T=T)
    umf4 <- unmarkedFramePCO(y = y.sim4, numPrimary=T)
    m4 <- pcountOpen(~1, ~1, ~1, ~1, umf4, K=30, dynamics="autoreg",
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m4)
    simout4[i, 1:2] <- exp(e[1:2])
    simout4[i, 3:4] <- plogis(e[3:4])
    }

png("pcountOpenSim4.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout4[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout4[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout4[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout4[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()










# No trend model




sim5 <- function(lambda=1, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    gamma <- (1-omega)*lambda
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    y[] <- rbinom(M*T, N, p)
    return(y)
}



set.seed(3223)
nsim5 <- 500
simout5 <- matrix(NA, nsim5, 3)
colnames(simout5) <- c('lambda', 'omega', 'p')
for(i in 1:nsim5) {
    cat("sim5", i, "\n"); flush.console()
    lambda <- 1
    omega <- 0.7
    p <- 0.7
    T <- 5
    y.sim5 <- sim5(lambda, omega, p, T=T)
    umf5 <- unmarkedFramePCO(y = y.sim5, numPrimary=T)
    m5 <- pcountOpen(~1, ~1, ~1, ~1, umf5, K=20, dynamics="notrend",
        starts=c(log(lambda), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m5)
    simout5[i, 1] <- exp(e[1])
    simout5[i, 2:3] <- plogis(e[2:3])
    }

png("pcountOpenSim5.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout5[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout5[,2], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout5[,3], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()












## Simulate data with some dates[i,1] > 1
set.seed(333)

sim6 <- function(lambda=4, gamma=0.1, omega=0.8, p=0.7, M=100, T=5)
{
    y <- N <- date <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    date[,1] <- pmax(rpois(M, 3), 1)

    for(i in 1:M) {
    if(date[i,1] > 1) {
        for(d in 2:date[i, 1]) {
            S1 <- rbinom(1, N[i, 1], omega)
            G1 <- rpois(1, gamma)
            N[i, 1] <- S1 + G1
            }
        }
    for(t in 2:T) {
        date[i, t] <- date[i, t-1] + 1
        S[i, t-1] <- rbinom(1, N[i, t-1], omega)
        G[i, t-1] <- rpois(1, gamma)
        N[i, t] <- S[i, t-1] + G[i, t-1]
        }}
    y[] <- rbinom(M*T, N, p)
    mode(date) <- "integer"
    return(list(y=y, dates=date))
}




set.seed(3223)
nsim6 <- 50
simout6 <- matrix(NA, nsim6, 4)
colnames(simout6) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim6) {
    cat("sim6", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    T <- 5
    yd <- sim6(lambda, gamma, omega, p, M=100, T=T)
    y.sim6 <- yd$y
    dates6 <- yd$dates
    umf6 <- unmarkedFramePCO(y = y.sim6, primaryPeriod=dates6, numPrimary=T)
    m6 <- pcountOpen(~1, ~1, ~1, ~1, umf6, K=15,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m6)
    simout6[i, 1:2] <- exp(e[1:2])
    simout6[i, 3:4] <- plogis(e[3:4])
    }



png("pcountOpenSim6.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout6[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout6[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout6[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout6[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()














## Simulate no covariates, constant sampling period intervals,
## WITH secondary samples

sim7 <- function(lambda=1, gamma=0.5, omega=0.8, p=0.7, M=100, T=5, J=3)
{
    y <- matrix(NA, M, J*T)
    N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega)
        G[,t] <- rpois(M, gamma)
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, p)
    return(y)
}




set.seed(3223)
nsim7 <- 500
simout7 <- matrix(NA, nsim7, 4)
colnames(simout7) <- c('lambda', 'gamma', 'omega', 'p')
for(i in 1:nsim7) {
    cat("sim7", i, "\n"); flush.console()
    lambda <- 1
    gamma <- 0.5
    omega <- 0.8
    p <- 0.7
    T <- 5
    y.sim7 <- sim7(lambda, gamma, omega, p, T=T)
    umf7 <- unmarkedFramePCO(y = y.sim7, numPrimary=T)
    m7 <- pcountOpen(~1, ~1, ~1, ~1, umf7, K=15,
        starts=c(log(lambda), log(gamma), plogis(omega), plogis(p)), se=FALSE)
    e <- coef(m7)
    simout7[i, 1:2] <- exp(e[1:2])
    simout7[i, 3:4] <- plogis(e[3:4])
    }

png("pcountOpenSim7.png", width=6, height=6, units="in", res=360)
par(mfrow=c(2,2))
hist(simout7[,1], xlab=expression(lambda)); abline(v=lambda, lwd=2, col=4)
hist(simout7[,2], xlab=expression(gamma)); abline(v=gamma, lwd=2, col=4)
hist(simout7[,3], xlab=expression(omega)); abline(v=omega, lwd=2, col=4)
hist(simout7[,4], xlab=expression(p)); abline(v=p, lwd=2, col=4)
dev.off()















## Simulate covariate model with constant intervals and secondary samples



sim8 <- function(lam=c(0,1), gam=c(-1,-1), om=c(2,-1), p=c(-1,1), M=100,
    T=5, J=3)
{
    y <- det <- matrix(NA, M, J*T)
    gamma <- omega <- N <- matrix(NA, M, T)
    S <- G <- matrix(NA, M, T-1)
    veght <- rnorm(M)
    isolation <- matrix(rnorm(M*T), M, T)
    time <- matrix(rnorm(M*J*T, 1), M, J*T)
    lambda <- exp(lam[1] + lam[2]*veght)
    gamma[] <- exp(gam[1] + gam[2]*isolation)
    omega[] <- plogis(om[1] + om[2]*isolation)
    det[] <- plogis(p[1] + p[2]*time)

    N[,1] <- rpois(M, lambda)
    for(t in 1:(T-1)) {
        S[,t] <- rbinom(M, N[,t], omega[,t])
        G[,t] <- rpois(M, gamma[,t])
        N[,t+1] <- S[,t] + G[,t]
        }
    N <- N[,rep(1:T, each=J)]
    y[] <- rbinom(M*J*T, N, det)
    return(list(y=y, covs=data.frame(veght=veght,
        isolation=isolation, time=time)))
}





nsim8 <- 500
simout8 <- matrix(NA, nsim8, 8)
colnames(simout8) <- c('lam0', 'lam1', 'gam0', 'gam1', 'om0', 'om1', 'p0', 'p1')
for(i in 1:nsim8) {
    cat("sim8", i, "\n"); flush.console()
    lam <- c(-2, 1)
    gam <- c(-1, -1)
    om <- c(0, -1)
    p <- c(-1, 1)
    T <- 5
    sim8out <- sim8(lam, gam, om, p, T=T)
    y.sim8 <- sim8out$y
    covs <- sim8out$covs
    cn <- colnames(covs)
    siteCovs <- covs[,grep("veght", cn), drop=FALSE]
    yearlySiteCovs <- list(isolation=covs[,grep("isolation", cn)])
    obsCovs <- list(time = covs[,grep("time", cn)])
    umf8 <- unmarkedFramePCO(y = y.sim8, siteCovs=siteCovs,
        yearlySiteCovs=yearlySiteCovs, obsCovs=obsCovs, numPrimary=T)
    m8 <- pcountOpen(~veght, ~isolation, ~isolation, ~time, umf8, K=30, se=F,
        starts=c(lam, gam, om, p))
    e <- coef(m8)
    simout8[i, ] <- e
    }

png("pcountOpenSim8.png", width=6, height=8, units="in", res=360)
par(mfrow=c(4,2))
hist(simout8[,1], xlab=expression(lambda)); abline(v=lam[1], lwd=2, col=4)
hist(simout8[,2], xlab=expression(lambda)); abline(v=lam[2], lwd=2, col=4)
hist(simout8[,3], xlab=expression(gamma)); abline(v=gam[1], lwd=2, col=4)
hist(simout8[,4], xlab=expression(gamma)); abline(v=gam[2], lwd=2, col=4)
hist(simout8[,5], xlab=expression(omega)); abline(v=om[1], lwd=2, col=4)
hist(simout8[,6], xlab=expression(omega)); abline(v=om[2], lwd=2, col=4)
hist(simout8[,7], xlab=expression(p)); abline(v=p[1], lwd=2, col=4)
hist(simout8[,8], xlab=expression(p)); abline(v=p[2], lwd=2, col=4)
dev.off()









