cooccu <- function(psiformulaA = ~1, psiformulaB = ~1, dependformula = ~1,
                   pformulaA = ~1, pformulaB = ~1,
                   fpformulaA = NULL, fpformulaB = NULL,
                   data,
                   interaction=c("dominance", "general"),
                   fpmodel=c("full", "confusion"),
                   fix=NULL, # Could be "gamma", ???
                   starts, method="BFGS", se=TRUE, ...) {

interaction <- match.arg(interaction)
fpmodel <- match.arg(fpmodel)

# Could loop over length(D) and assign() each component to object
yA <- D$yA
yB <- D$yB
XpsiA <- D$XpsiA
XpsiB <- D$XpsiB
Xgamma <- D$Xgamma
XpA <- D$XpA
XpB <- D$XpB
XfpA <- D$XfpA
XfpB <- D$XfpB
XpsiA.offset <- D$XpsiA.offset # getDesign should return 0's if no offset
XpsiB.offset <- D$XpsiB.offset
Xgamma.offset <- D$Xgamma.offset
XpA.offset <- D$XpA.offset
XpB.offset <- D$XpB.offset
XfpA.offset <- D$XfpA.offset
XfpB.offset <- D$XfpB.offset
removed <- D$removed.sites
noFPsites <- D$noFPsites # Should be logical. TRUE if Pr(fp)=0

R <- nrow(y)
J <- ncol(y)
nSpp <- 2

nOPA <- ncol(XpsiA)
nOPB <- ncol(XpsiB)
nGP <- ncol(Xgamma)
nPA <- ncol(XpA)
nPB <- ncol(XpB)
nFPA <- ncol(XfpA)
nFPB <- ncol(XfpB)
nP <- nOPA + nOPB + nGP + nPA + nPB + nFPA + nFPB
if(missing(starts))
    starts <- rep(0, nP)
else if(length(starts) != nP)
    stop("There should be ", nP, "starting values, not", length(starts))


opaNames <- colnames(XpsiA)
opbNames <- colnames(XpsiB)
gpNames <- colnames(Xgamma)
paNames <- colnames(XpA)
pbNames <- colnames(XpB)
fpaNames <- colnames(XfpA)
fpbNames <- colnames(XfpB)


nll <- function(pars) {
    psiA <- plogis(XpsiA %*% pars[1:nOPA] + XpsiA.offset)
    psiB <- plogis(XpsiB %*% pars[(nOPA+1):(nOPA+nOPB)] + XpsiB.offset)
    if(identical(interaction, "dominance")) {
        gamma <- plogis(Xgamma %*% pars[(nOPA+nOPB+1):(nGP)] +
                        Xgamma.offset)
        psiBz <- psiB*gamma
        psiAB <- psiA*psiBz
        phi <- cbind(psiAB, psiA*(1-psiBz), psiB*(1-psiA),
                     (1-psiA)*(1-psiB))
    } else if(identical(interaction, "general")) {
        gamma <-  exp(Xgamma %*% pars[(nOPA+nOPB+1):(nGP)] + Xgamma.offset)
        psiAB <- psiA * psiB * gamma
        if(any(pmax(psiA + psiB - 1, 0) > psiAB) |
           any(psiAB > pmin(psiA, psiB))) {
            # warning("Model constraints violated.")
            return(NA)
        }
        phi <- cbind(psiAB, psiA-psiAB, psiB-psiAB, 1-psiA-psiB-psiAB)
    }
    if(any(rowSums(phi) < 0.99))
        stop("sum(phi) less than 0.99.")
    # Detection prob of A, possibly dependent upon presence of B
    pApars <- pApars.B0 <- pars[(nOPA+nOPB+nGP+1):(nOPA+nOPB+nGP+nPA)]
    pA.B1 <- plogis(XpA %*% pApars + XpA.offset)
    # Do the next 2 lines if SppB is in pformulaA, otherwise pA.B0=pA.B1
    pApars.B0["SppB"] <- 0
    pA.B0 <- plogis(XpA %*% pApars.B0 + XpA.offset)
    # Detection prob of B, possibly dependent upon presence of A
    pBpars <- pBpars.A0 <- pars[(nOPA+nOPB+nGP+nPA+1):
                                (nOPA+nOPB+nGP+nPA+nPB)]
    pBpars.A0["SppA"] <- 0
    pB.A1 <- plogis(XpB %*% pBpars + XpB.offset)
    pB.A0 <- plogis(XpB %*% pBpars.A0 + XpB.offset)
    # Only do the next 2 steps if fpformulae are not NULL
    # False-positive prob of A, possibly dependent upon presence of B
    fpApars <- fpApars.B0 <- pars[(nOPA+nOPB+nGP+nPA+nPB+1):
                                  (nOPA+nOPB+nGP+nPA+nPB+nFPA)]
    fpApars.B0["SppB"] <- 0
    fpB.A1 <- plogis(XfpA %*% fpApars + XfpA.offset)
    fpB.A0 <- plogis(XfpA %*% fpApars.B0 + XfpA.offset)
    # False-positive prob of B, possibly dependent upon presence of A
    fpBpars <- fpBpars.A0 <- pars[(nOPA+nOPB+nGP+nPA+nPB+nFPA+1):
                                  (nOPA+nOPB+nGP+nPA+nPB+nFPA+nFPB)]
    fpBpars.A0["SppA"] <- 0
    fpB.A1 <- plogis(XfpB %*% fpBpars + XfpB.offset)
    fpB.A0 <- plogis(XfpB %*% fpBpars.A0 + XfpB.offset)
    # Conditional on z, observation probs
    bin.A1.B1 <- dbinom(yA, 1, pA.B1, log=TRUE) +
        dbinom(yB, 1, pB.A1, log=TRUE)
    bin.A1.B0 <- dbinom(yA, 1, pA.B0, log=TRUE) +
        dbinom(yB, 1, fpB.A0, log=TRUE)
    bin.A0.B1 <- dbinom(yA, 1, fpA.B1, log=TRUE) +
        dbinom(yB, 1, pB.A1, log=TRUE)
    if(identical(fpmodel, "full"))
        bin.A0.B0 <- dbinom(yA, 1, fpA.B0, log=TRUE) +
            dbinom(yB, 1, fpB.A0, log=TRUE)
    else if(identical(fpmodel, "confusion"))
        bin.A0.B0 <- dbinom(yA, 1, 0, log=TRUE) +
            dbinom(yB, 1, 0, log=TRUE)
    # Handle NAs here
    mu <- cbind(exp(rowSums(bin.A1.B1)),
                exp(rowSums(bin.A1.B0)),
                exp(rowSums(bin.A0.B1)),
                exp(rowSums(bin.A0.B0)))
    L <- rowSums(phi * mu)
    -sum(log(L))
   }

fm <- optim(starts, nll, method=method, hessian=se, ...)
opt <- fm
if(se)
    covMat <- solve(fm$hessian) # Add try() or tryCatch()
else
    covMat <- matrix(NA, nP, nP)
mle <- fm$par
fmAIC <- 2 * fm$value + 2 * nP
names(mle) <- c(opaNames, opbNames, gpNames, paNames, pbNames,
                fpaNames, fpbNames)

psiA <- unmarkedEstimate(name = "Occupancy species A", short.name="psiA",
                         estimates = mle[1:nOPA],
                         covMat = covMat[1:nOPA, 1:nOPA, drop=FALSE],
                         invlink = "logistic",
                         invlinkGrad = "logistic.grad")
psiB <- unmarkedEstimate(name = "Occupancy species B", short.name="psiB",
                         estimates = mle[(nOPA+1):(nOPA+nOPB)],
                         covMat = covMat[(nOPA+1):(nOPA+nOPB),
                                         (nOPA+1):(nOPA+nOPB), drop=FALSE],
                         invlink = "logistic",
                         invlinkGrad = "logistic.grad")
gamma <- unmarkedEstimate(name = "Interaction", short.name="psiB",
                         estimates = mle[(nOPA+nOPB+1):(nOPA+nOPB+nGP)],
                         covMat = covMat[(nOPA+nOPB+1):(nOPA+nOPB+nGP),
                                         (nOPA+nOPB+1):(nOPA+nOPB+nGP),
                                         drop=FALSE],
                         invlink = if(interaction=="dominance")
                                      "logistic"
                                   else "exp",
                         invlinkGrad = if(interaction=="dominance")
                                          "logistic.grad"
                                       else "exp")
pA <- unmarkedEstimate(name = "False negative species A", short.name="pA",
                       estimates = mle[(nOPA+nOPB+nGP+1):
                                       (nOPA+nOPB+nGP+nPA)],
                       covMat = covMat[(nOPA+nOPB+nGP+1):
                                       (nOPA+nOPB+nGP+nPA),
                                       (nOPA+nOPB+nGP+1):
                                       (nOPA+nOPB+nGP+nPA), drop=FALSE],
                       invlink = "logistic", invlinkGrad = "logistic.grad")
pB <- unmarkedEstimate(name = "False negative species B", short.name="pB",
                       estimates = mle[(nOPA+nOPB+nGP+nPA+1):
                                       (nOPA+nOPB+nGP+nPA+nPB)],
                       covMat = covMat[(nOPA+nOPB+nGP+nPA+1):
                                       (nOPA+nOPB+nGP+nPA+nPB),
                                       (nOPA+nOPB+nGP+nPA+1):
                                       (nOPA+nOPB+nGP+nPA+nPB),
                                       drop=FALSE],
                       invlink = "logistic", invlinkGrad = "logistic.grad")
estimateList <- unmarkedEstimateList(list(psiA=psi, psiB=psiB,
                                          gamma=gamma, pA=pA, pB=pB))
if(fpA) {
    estimateList@estimates$fpA <-
        unmarkedEstimate(name = "False positive species A",
                         short.name="pA",
                         estimates = mle[(nOPA+nOPB+nGP+nPA+nPB+1):
                                         (nOPA+nOPB+nGP+nPA+nPB+nFPA)],
                       covMat = covMat[(nOPA+nOPB+nGP+nPA+nPB+1):
                                       (nOPA+nOPB+nGP+nPA+nPB+nFPA),
                                       (nOPA+nOPB+nGP+nPA+nPB+1):
                                       (nOPA+nOPB+nGP+nPA+nPB+nFPA),
                                       drop=FALSE],
                       invlink = "logistic", invlinkGrad = "logistic.grad")
}
if(fpB) {
    estimateList@estimates$fpB <-
        unmarkedEstimate(name = "False positive species B",
                         short.name="pB",
                         estimates = mle[(nOPA+nOPB+nGP+nPA+nPB+nFPA+1):
                                        (nOPA+nOPB+nGP+nPA+nPB+nFPA+nFPB)],
                       covMat = covMat[(nOPA+nOPB+nGP+nPA+nPB+nFPA+1):
                                       (nOPA+nOPB+nGP+nPA+nPB+nFPA+nFPB),
                                       (nOPA+nOPB+nGP+nPA+nPB+nFPA+1):
                                       (nOPA+nOPB+nGP+nPA+nPB+nFPA+nFPB),
                                       drop=FALSE],
                       invlink = "logistic", invlinkGrad = "logistic.grad")
}

umfit <- new("unmarkedFitCooccu", fitType=""
             call=match.call(),




}

