


e2dist1 <- function (x, y)
{
    i <- sort(rep(1:nrow(y), nrow(x)))
    dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
    matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}








sim.data <- function(N=25, T=5, sigma=0.5, lam0=1, knownID=0, lc=10,
    buff=3, sslower=0, ssupper=15)
{
    coords <- seq(sslower+buff, ssupper-buff, length=lc)
    X <- cbind(x=rep(coords, each=lc), y=rep(coords, times=lc))
    npts <- nrow(X)
    # Home range centers
    sx <- runif(N, sslower, ssupper)
    sy <- runif(N, sslower, ssupper)
    S <- cbind(sx, sy)
    D <- e2dist1(S, X)
    lam <- lam0*exp(-(D*D)/(2*sigma*sigma))
    Z <- array(NA, c(N, npts, T))
    for(t in 1:T) {
        Z[,,t] <- rpois(prod(dim(D)), lam)
    }
    y <- apply(Z, c(2,3), sum)
    capind <- (1:N)[apply(Z>0, 1, any)]
    if(length(capind) < knownID)
        stop("try again. Some known guys were never caught")
    Zknown <- NULL
    if(knownID>0) {
        known.kp <- capind[1:knownID]
        Zknown <- Z[known.kp,,]
    }

    list(y=y, N=N, T=T, sigma=sigma, lam0=lam0, X=X, S=S, buff=buff,
         sslower=sslower, ssupper=ssupper, knownID=knownID,
         Z=Z, Zknown=Zknown)

}






funcZ <- function(y, X, Zknown, M, niters, sslower=0, ssupper=15) {

    R <- nrow(y)
    T <- ncol(y)
    S <- cbind(runif(M, sslower, ssupper), runif(M, sslower, ssupper))
    D <- e2dist1(S, X)
    theta <- runif(1, .25, 1.25)
    lam0 <- runif(1, .1, 1)
    lam <- lam0*exp(-(D*D)/(2*theta*theta))
    w <- rbinom(M,1,.5)
    psi <- runif(1,.2,.8)

    Z <- array(NA, c(M,R,T))
    nMarked <- 0
    marked <- rep(FALSE, M)
    if(!missing(Zknown) || !is.null(Zknown)) {
        nMarked <- nrow(Zknown)
        marked[1:nMarked] <- TRUE
        Z[1:nMarked,,] <- Zknown
    }
    w[marked] <- 1
    Zdata <- !is.na(Z)
    for(r in 1:R) {
        for(t in 1:T) {
            if(y[r,t]==0) {
                Z[,r,t] <- 0
                next
            }
            unmarked <- !Zdata[,r,t]
            nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
            if(nUnknown < 0)
                browser()
            probs <- lam[,r]*w
            probs <- probs[unmarked]
            Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs)
        }
    }

    out <- matrix(NA,nrow=niters,ncol=4)
    colnames(out) <- c("sigma", "lam0", "psi", "N")

    cat("\nstarting values =", c(theta, lam0, psi, sum(w)), "\n\n")

    for(iter in 1:niters) {

        if(iter %% 100 == 0) {
            cat("iter", iter, format(Sys.time(), "%H:%M:%S"), "\n")
            cat("   current =", out[iter-1,], "\n")
        }

        # update theta
        theta.cand <- rnorm(1, theta, 0.03)
        if(theta.cand>0){
            lam.cand <- lam0*exp(-(D*D)/(2*theta.cand*theta.cand))
            # w is recycled over lam, R times
            # lam*w is recycled over Z, T times
            ll <- sum(dpois(Z, lam*w, log=TRUE))
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            if(runif(1) < exp( llcand  - ll ) ){
                ll<-llcand
                lam<-lam.cand
                theta<-theta.cand
            }
        }

        # update lam0
        lam0.cand <- rnorm(1, lam0, .06)
        if(lam0.cand>0) {
            lam.cand <- lam0.cand*exp(-(D*D)/(2*theta*theta))
            llcand <- sum(dpois(Z, lam.cand*w, log=TRUE))
            if(runif(1) < exp(llcand - ll)) {
                ll <- llcand
                lam0 <- lam0.cand
                lam <- lam.cand
            }
        }

        ### update "w" here
        wUps <- 0
        seen <- apply(Z>0, 1, any)
        for(i in 1:M) {
            if(seen[i] | marked[i])
                next
            wcand <- ifelse(w[i]==0, 1, 0)

            ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
            llcand <- sum(dpois(Z[i,,], lam[i,]*wcand, log=TRUE))

            prior <- dbinom(w[i], 1, psi, log=TRUE)
            prior.cand <- dbinom(wcand, 1, psi, log=TRUE)
            if(runif(1) < exp( (llcand+prior.cand) - (ll+prior) )) {
                w[i] <- wcand
                wUps <- wUps+1
            }
        }

        # update Z
        for(r in 1:R) {
            zip <- lam[,r]*w
            probs <- zip/sum(zip)
            for(t in 1:T) {
                if(y[r,t]==0) {
                    Z[,r,t] <- 0
                    next
                }
                unmarked <- !Zdata[,r,t]
                nUnknown <- y[r,t] - sum(Z[!unmarked,r,t])
                if(nUnknown < 0)
                    browser()
                Z[unmarked,r,t] <- rmultinom(1, nUnknown, probs[unmarked])
            }
        }


        # update psi
        psi<-rbeta(1,1+sum(w),1+M-sum(w))

        # update S
        Sups <- 0
        for(i in 1:M) {   # note this is "M" in general
            Scand <- c(rnorm(1, S[i,1], 0.1),
                       rnorm(1, S[i,2], 0.1))
            inbox <- Scand[1]>=sslower & Scand[1]<=ssupper &
                     Scand[2]>=sslower & Scand[2]<=ssupper
            if(inbox) {
                dtmp <- sqrt((Scand[1] - X[,1])^2 + (Scand[2] - X[,2])^2)
                lam.cand <- lam0*exp(-(dtmp*dtmp)/(2*theta*theta) )

                # recycle lam over Z
                ll <- sum(dpois(Z[i,,], lam[i,]*w[i], log=TRUE))
                llcand <- sum(dpois(Z[i,,], lam.cand*w[i], log=TRUE))

                if(runif(1) < exp(llcand - ll)) {
                    ll <- llcand
                    S[i,] <- Scand
                    lam[i,] <- lam.cand
                    D[i,] <- dtmp
                    Sups <- Sups+1
                }
            }
        }

        if(iter %% 100 == 0) {
            cat("   Acceptance rates\n")
            cat("     w =", wUps/M, "\n")
            cat("     S =", Sups/M, "\n")
        }

        out[iter,] <- c(theta,lam0,psi,sum(w) )

    }

    return(out)
}









# Example
# library(coda)
# str(sim1 <- sim.data(N=25, T=5, sigma=1, lam0=1, knownID=10))
# fm1 <- funcZ(sim1$y, sim1$X, sim1$Zknown, M=50, niters=6000)
# plot(mcmc(fm1))
# rejectionRate(window(mcmc(fm1), start=3001))
