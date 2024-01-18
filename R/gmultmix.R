
# data will need to be an unmarkedMultFrame
gmultmix <- function(lambdaformula, phiformula, pformula, data,
    mixture=c("P", "NB", "ZIP"), K, starts, method = "BFGS", se = TRUE,
    engine=c("C","R"), threads=1, ...)
{
if(!is(data, "unmarkedFrameGMM"))
    stop("Data is not of class unmarkedFrameGMM.")

engine <- match.arg(engine, c("C", "R"))

mixture <- match.arg(mixture)

formlist <- list(lambdaformula = lambdaformula, phiformula = phiformula,
    pformula = pformula)
check_no_support(formlist)
form <- as.formula(paste(unlist(formlist), collapse=" "))
D <- getDesign(data, formula = form)

Xlam <- D$Xlam
Xphi <- D$Xphi
Xdet <- D$Xdet
y <- D$y  # MxJT

Xlam.offset <- D$Xlam.offset
Xphi.offset <- D$Xphi.offset
Xdet.offset <- D$Xdet.offset
if(is.null(Xlam.offset)) Xlam.offset <- rep(0, nrow(Xlam))
if(is.null(Xphi.offset)) Xphi.offset <- rep(0, nrow(Xphi))
if(is.null(Xdet.offset)) Xdet.offset <- rep(0, nrow(Xdet))

if(missing(K) || is.null(K)) K <- max(y, na.rm=TRUE) + 100
k <- 0:K
lk <- length(k)
M <- nrow(y)
T <- data@numPrimary
R <- numY(data) / T
J <- obsNum(data) / T

y <- array(y, c(M, R, T))
y <- aperm(y, c(1,3,2))
yt <- apply(y, 1:2, function(x) {
    if(all(is.na(x)))
        return(NA)
    else return(sum(x, na.rm=TRUE))
    })


piFun <- data@piFun

lamPars <- colnames(Xlam)
detPars <- colnames(Xdet)
nLP <- ncol(Xlam)
if(T==1) {
    nPP <- 0
    phiPars <- character(0)
} else if(T>1) {
    nPP <- ncol(Xphi)
    phiPars <- colnames(Xphi)
    }
nDP <- ncol(Xdet)
nP <- nLP + nPP + nDP + (mixture%in%c('NB','ZIP'))
if(!missing(starts) && length(starts) != nP)
    stop(paste("The number of starting values should be", nP))


lfac.k <- lgamma(k+1)
kmyt <- array(NA, c(M, T, lk))
lfac.kmyt <- array(0, c(M, T, lk))
fin <- matrix(NA, M, lk)
naflag <- array(NA, c(M, T, R))
for(i in 1:M) {
    fin[i, ] <- k - max(yt[i,], na.rm=TRUE) >= 0
    for(t in 1:T) {
        naflag[i,t,] <- is.na(y[i,t,])
        if(!all(naflag[i,t,])) {
            kmyt[i,t,] <- k - yt[i,t]
            lfac.kmyt[i, t, fin[i,]] <- lgamma(kmyt[i, t, fin[i,]] + 1)
            }
        }
    }

nll_R <- function(pars) {
    lambda <- exp(Xlam %*% pars[1:nLP] + Xlam.offset)
    if(T==1)
        phi <- 1
    else if(T>1)
        phi <- drop(plogis(Xphi %*% pars[(nLP+1):(nLP+nPP)] + Xphi.offset))
    p <- plogis(Xdet %*% pars[(nLP+nPP+1):(nLP+nPP+nDP)] + Xdet.offset)

    phi.mat <- matrix(phi, M, T, byrow=TRUE)
    phi <- as.numeric(phi.mat)

    p <- matrix(p, nrow=M, byrow=TRUE)
    p <- array(p, c(M, J, T))
    p <- aperm(p, c(1,3,2))
    cp <- array(as.numeric(NA), c(M, T, R+1))

    for(t in 1:T) cp[,t,1:R] <- do.call(piFun, list(p[,t,]))
    cp[,,1:R] <- cp[,,1:R] * phi
    cp[,, 1:R][is.na(y)]<- NA   # andy added 5/29
    cp[,,R+1] <- 1 - apply(cp[,,1:R,drop=FALSE], 1:2, sum, na.rm=TRUE)

    switch(mixture,
      P = f <- sapply(k, function(x) dpois(x, lambda)),
      NB = f <- sapply(k, function(x) dnbinom(x, mu=lambda, size=exp(pars[nP]))),
      ZIP = f <- sapply(k, function(x) dzip(rep(x, length(lambda)), lambda=lambda, psi=plogis(pars[nP])))
    )
    g <- matrix(as.numeric(NA), M, lk)
    for(i in 1:M) {
        A <- matrix(0, lk, T)
        for(t in 1:T) {
            na <- naflag[i,t,]
            if(!all(na))
                A[, t] <- lfac.k - lfac.kmyt[i, t,] +
                    sum(y[i, t, !na] * log(cp[i, t, which(!na)])) +
                    kmyt[i, t,] * log(cp[i, t, R+1])
            }
        g[i,] <- exp(rowSums(A))
        }
    f[!fin] <- g[!fin] <- 0
    ll <- rowSums(f*g)
    -sum(log(ll))
    }


if(engine=="R"){
  nll <- nll_R
} else {
  long_format <- function(x){
    out <- matrix(aperm(x,c(1,3,2)),nrow=nrow(x),ncol=dim(x)[2]*dim(x)[3])
    as.vector(t(out))
  }
  y_long <- long_format(y)
  kmytC <- kmyt
  kmytC[which(is.na(kmyt))] <- 0

  mixture_code <- switch(mixture, P={1}, NB={2}, ZIP={3})
  n_param <- c(nLP, nPP, nDP, mixture%in%c("NB","ZIP"))
  Kmin <- apply(yt, 1, max, na.rm=TRUE)

  nll <- function(params) {
    nll_gmultmix(params, n_param, y_long, mixture_code, piFun, Xlam, Xlam.offset,
                 Xphi, Xphi.offset, Xdet, Xdet.offset, k, lfac.k, lfac.kmyt,
                 kmytC, Kmin, threads)
  }

  if(!piFun%in%c('doublePiFun','removalPiFun','depDoublePiFun')){
    warning("Custom pi functions are not supported by C engine. Using R engine instead.")
    nll <- nll_R
  }

}

if(missing(starts)) starts <- rep(0, nP)
fm <- optim(starts, nll, method = method, hessian = se, ...)
covMat <- invertHessian(fm, nP, se)
ests <- fm$par
fmAIC <- 2 * fm$value + 2 * nP

nbParm <- switch(mixture, P={character(0)}, NB={"alpha"}, ZIP={"psi"})

names(ests) <- c(lamPars, phiPars, detPars, nbParm)

lamEstimates <- unmarkedEstimate(name = "Abundance", short.name = "lambda",
    estimates = ests[1:nLP],
    covMat = as.matrix(covMat[1:nLP, 1:nLP]), invlink = "exp",
    invlinkGrad = "exp")
estimateList <- unmarkedEstimateList(list(lambda=lamEstimates))

if(T>1) {
    phiEstimates <- unmarkedEstimate(name = "Availability",
                                     short.name = "phi",
                                     estimates = ests[(nLP+1):(nLP+nPP)],
                                     covMat = as.matrix(covMat[(nLP+1) :
                                       (nLP+nPP), (nLP+1):(nLP+nPP)]),
                                     invlink = "logistic",
                                     invlinkGrad = "logistic.grad")
    estimateList@estimates$phi <- phiEstimates
}

detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
    estimates = ests[(nLP+nPP+1):(nLP+nPP+nDP)],
    covMat = as.matrix(
        covMat[(nLP+nPP+1):(nLP+nPP+nDP), (nLP+nPP+1):(nLP+nPP+nDP)]),
    invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList@estimates$det <- detEstimates

if(identical(mixture,"NB"))
    estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
        short.name = "alpha", estimates = ests[nP],
        covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
        invlinkGrad = "exp")

if(identical(mixture,"ZIP")) {
    estimateList@estimates$psi <- unmarkedEstimate(name="Zero-inflation",
        short.name = "psi", estimates = ests[nP],
        covMat=as.matrix(covMat[nP, nP]), invlink = "logistic",
        invlinkGrad = "logistic.grad")
}

umfit <- new("unmarkedFitGMM", fitType = "gmn",
    call = match.call(), formula = form, formlist = formlist,
    data = data, estimates = estimateList, sitesRemoved = D$removed.sites,
    AIC = fmAIC, opt = fm, negLogLike = fm$value, nllFun = nll,
    mixture=mixture, K=K)

return(umfit)
}




