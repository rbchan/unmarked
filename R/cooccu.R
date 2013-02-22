cooccu <- function(psiformulaA = ~1, psiformulaB = ~1,
                   gammaformula = NULL,
                   pformulaA = ~1, pformulaB = ~1,
                   fpformulaA = NULL, fpformulaB = NULL,
                   data,
                   interaction=c("dominance", "general"),
                   fpmodel=c("full", "confusion"),
                   starts, method="BFGS", se=TRUE, ...) {

interaction <- match.arg(interaction)
fpmodel <- match.arg(fpmodel)

if(identical(fpmodel, "confusion")) {
    if("SppB" %in% all.vars(fpformulaA))
        stop("It doesn't make sense to have 'SppB' in fpformulaA when fpmodel='confusion'")
    if("SppA" %in% all.vars(fpformulaB))
        stop("It doesn't make sense to have 'SppA' in fpformulaB when fpmodel='confusion'")
}

formlist <- list(psiformulaA=psiformulaA, psiformulaB=psiformulaB,
                 gammaformula=gammaformula,
                 pformulaA=pformulaA, pformulaB=pformulaB,
                 fpformulaA=fpformulaA, fpformulaB=fpformulaB)
formula <- as.formula(paste(unlist(formlist), collapse=" "))

#D <- getDesign(data, formula)
D <- getDesign(data, psiformulaA, psiformulaB, gammaformula,
               pformulaA, pformulaB, fpformulaA, fpformulaB)
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
noFP <- D$noFP # Should be logical. TRUE if Pr(fp)=0

if(is.null(psiformulaA) | is.null(psiformulaB) |
   is.null(pformulaA) | is.null(pformulaB))
    stop("Only gammaformulaA, fpformulaA, and fpformulaB can be NULL")
fixgamma <- if(is.null(gammaformula)) TRUE else FALSE
fixfpA <- if(is.null(fpformulaA)) TRUE else FALSE
fixfpB <- if(is.null(fpformulaB)) TRUE else FALSE

R <- nrow(yA)
J <- ncol(yA)
nSpp <- 2

nOPA <- ncol(XpsiA)
nOPB <- ncol(XpsiB)
nGP <- ncol(Xgamma)
nPA <- ncol(XpA)
nPB <- ncol(XpB)
nFPA <- ncol(XfpA)
nFPB <- ncol(XfpB)
nP <- nOPA + nOPB + nGP + nPA + nPB + nFPA + nFPB

opaNames <- colnames(XpsiA)
opbNames <- colnames(XpsiB)
gpNames <- colnames(Xgamma)
paNames <- colnames(XpA)
pbNames <- colnames(XpB)
fpaNames <- colnames(XfpA)
fpbNames <- colnames(XfpB)

if("SppA" %in% c(opaNames, opbNames, gpNames, paNames, fpaNames))
    stop("The variable 'SppA' can only be in pformulaB or fpformulaB")
if("SppB" %in% c(opaNames, opbNames, gpNames, pbNames, fpbNames))
    stop("The variable 'SppB' can only be in pformulaA or fpformulaA")

if(missing(starts)) {
    starts <- rep(0, nP)
    # Start the false-positive probs (intercepts) at small value
    starts[nOPA+nOPB+nGP+nPA+nPB+1] <- -3
    starts[nOPA+nOPB+nGP+nPA+nPB+nFPA+1] <- -3
}
else if(length(starts) != nP)
    stop("There should be ", nP, " starting values, not ", length(starts))
names(starts) <- c(opaNames, opbNames, gpNames, paNames, pbNames,
                  fpaNames, fpbNames)

nll <- function(pars) {
    psiA <- plogis(XpsiA %*% pars[1:nOPA] + XpsiA.offset)
    psiB <- plogis(XpsiB %*% pars[(nOPA+1):(nOPA+nOPB)] + XpsiB.offset)
    if(identical(interaction, "dominance")) {
        if(fixgamma)
            gamma <- rep(1, R)
        else
            gamma <- plogis(Xgamma %*%pars[(nOPA+nOPB+1):(nOPA+nOPB+nGP)] +
                            Xgamma.offset)
        psiBz <- psiB*gamma
        psiAB <- psiA*psiBz
        phi <- cbind(psiAB, psiA*(1-psiBz), psiB*(1-psiA),
                     (1-psiA)*(1-psiB))
    } else if(identical(interaction, "general")) {
        if(fixgamma)
            gamma <- rep(1, R)
        else
            gamma <-  exp(Xgamma %*% pars[(nOPA+nOPB+1):(nOPA+nOPB+nGP)] +
                          Xgamma.offset)
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
    pA.B1 <- matrix(plogis(XpA %*% pApars + XpA.offset),
                    R, J, byrow=TRUE)
    if("SppB" %in% names(pApars)) {
        pApars.B0["SppB"] <- 0
        pA.B0 <- matrix(plogis(XpA %*% pApars.B0 + XpA.offset),
                        R, J, byrow=TRUE)
    } else pA.B0 <- pA.B1
    # Detection prob of B, possibly dependent upon presence of A
    pBpars <- pBpars.A0 <- pars[(nOPA+nOPB+nGP+nPA+1):
                                (nOPA+nOPB+nGP+nPA+nPB)]
    pB.A1 <- matrix(plogis(XpB %*% pBpars + XpB.offset),
                    R, J, byrow=TRUE)
    if("SppA" %in% names(pBpars)) {
        pBpars.A0["SppA"] <- 0
        pB.A0 <- matrix(plogis(XpB %*% pBpars.A0 + XpB.offset),
                        R, J, byrow=TRUE)
    } else pB.A0 <- pB.A1
    # False-positive prob of A, possibly dependent upon presence of B
    if(fixfpA) {
        fpA.B1 <- fpA.B0 <- matrix(0, R, J)
    } else {
#        browser()
        fpApars <- fpApars.B0 <- pars[(nOPA+nOPB+nGP+nPA+nPB+1):
                                      (nOPA+nOPB+nGP+nPA+nPB+nFPA)]
        fpA.B1 <- matrix(plogis(XfpA %*% fpApars + XfpA.offset),
                         R, J, byrow=TRUE)
        if("SppB" %in% names(fpApars)) {
            fpApars.B0["SppB"] <- 0
            fpA.B0 <- matrix(plogis(XfpA %*% fpApars.B0 + XfpA.offset),
                             R, J, byrow=TRUE)
        } else fpA.B0 <- fpA.B1
        if(any(noFP)) {
            fpA.B1[noFP] <- 0
            fpA.B0[noFP] <- 0
        }
    }
    # False-positive prob of B, possibly dependent upon presence of A
    if(fixfpB) {
        fpB.A1 <- fpB.A0 <- matrix(0, R, J)
    } else {
        fpBpars <- fpBpars.A0 <- pars[(nOPA+nOPB+nGP+nPA+nPB+nFPA+1):
                                      (nOPA+nOPB+nGP+nPA+nPB+nFPA+nFPB)]
        fpB.A1 <- matrix(plogis(XfpB %*% fpBpars + XfpB.offset),
                         R, J, byrow=TRUE)
        if("SppA" %in% names(fpBpars.A0)) {
            fpBpars.A0["SppA"] <- 0
            fpB.A0 <- matrix(plogis(XfpB %*% fpBpars.A0 + XfpB.offset),
                             R, J, byrow=TRUE)
        } else fpB.A0 <- fpB.A1
        if(any(noFP)) {
            fpB.A1[noFP] <- 0
            fpB.A0[noFP] <- 0
        }
    }
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
    # Handle NAs
    # Shouldn't be NAs in phi since those sites should've been removed
    bin.A1.B1[is.na(bin.A1.B1)] <- 0
    bin.A1.B0[is.na(bin.A1.B0)] <- 0
    bin.A0.B1[is.na(bin.A0.B1)] <- 0
    bin.A0.B0[is.na(bin.A0.B0)] <- 0
    mu <- cbind(exp(rowSums(bin.A1.B1)),
                exp(rowSums(bin.A1.B0)),
                exp(rowSums(bin.A0.B1)),
                exp(rowSums(bin.A0.B0)))
    L <- rowSums(phi * mu)
    -sum(log(L))
   }

fm <- optim(starts, nll, method=method, hessian=se, ...)
opt <- fm
if(se) {
    covMat <- try(solve(fm$hessian), silent=TRUE)
    if(identical(class(covMat), "try-error")) {
        warning("Unable to compute variance-covariance matrix because:\n",
                attr(covMat, "condition"),
                "Try providing good starting values or using fewer covariates",
                call.=FALSE)
        covMat <- matrix(NA, nP, nP)
    }
} else
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
estimateList <- unmarkedEstimateList(list(psiA=psiA, psiB=psiB,
                                          pA=pA, pB=pB))

if(!fixgamma) {
    estimateList@estimates$gamma <-
        unmarkedEstimate(name = "Interaction", short.name="psiB",
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
}


if(!fixfpA) {
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

if(!fixfpB) {
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

if(!fixgamma) {
    nms <- names(estimateList@estimates)
    nms <- nms[!nms %in% c("psiA", "psiB", "gamma")]
    estimateList@estimates <- estimateList@estimates[c("psiA", "psiB", "gamma", nms)]
}

umfit <- new("unmarkedFitCo", fitType="cooccu",
             call=match.call(), formula=formula,
             data=data, sitesRemoved = removed,
             estimates = estimateList, AIC = fmAIC, opt = opt,
             negLogLike = fm$value, nllFun = nll,
             interaction=interaction, fpmodel=fpmodel)

}

