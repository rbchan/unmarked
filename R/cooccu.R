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

yA <- D$yA
yB <- D$yB
XpsiA <- D$XpsiA
XpsiB <- D$XpsiB
Xgamma <- D$Xgamma
XpA <- D$XpA
XpB <- D$XpB
XfpA <- D$XfpA
XfpB <- D$XfpB


R <- nrow(y)
J <- ncol(y)



nll <- function(pars) {
    psiA <- plogis(XpsiA %*% pars[1:nOPA] + XpsiA.offset)
    psiB <- plogis(XpsiB %*% pars[(nOPA+1):(nOPA+nOPB)] + XpsiB.offset)
    if(identical(interaction, "dominance")) {
        gamma <- plogis(Xgamma %*% pars[(nOPA+nOPB+1):(nG)] +
                        Xgamma.offset)
        psiBz <- psiB*gamma
        psiAB <- psiA*psiBz
        phi <- cbind(psiAB, psiA*(1-psiBz), psiB*(1-psiA),
                     (1-psiA)*(1-psiB))
    } else if(identical(interaction, "general")) {
        gamma <-  exp(Xgamma %*% pars[(nOPA+nOPB+1):(nG)] + Xgamma.offset)
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
    pApars <- pApars.B0 <- pars[(nOPA+nOPB+nG+1):(nOPA+nOPB+nG+nPA)]
    pApars.B0["SppB"] <- 0
    pA.B1 <- plogis(XpA %*% pApars + XpA.offset)
    pA.B0 <- plogis(XpA %*% pApars.B0 + XpA.offset)
    # Detection prob of B, possibly dependent upon presence of A
    pBpars <- pBpars.A0 <- pars[(nOPA+nOPB+nG+nPA+1):
                                (nOPA+nOPB+nG+nPA+nPB)]
    pBpars.A0["SppA"] <- 0
    pB.A1 <- plogis(XpB %*% pBpars + XpB.offset)
    pB.A0 <- plogis(XpB %*% pBpars.A0 + XpB.offset)
    # Only do the next 2 steps if fpformulae are not NULL
    # False-positive prob of A, possibly dependent upon presence of B
    fpApars <- fpApars.B0 <- pars[(nOPA+nOPB+nG+nPA+nPB+1):
                                  (nOPA+nOPB+nG+nPA+nPB+nFPA)]
    fpApars.B0["SppB"] <- 0
    fpB.A1 <- plogis(XfpA %*% fpApars + XfpA.offset)
    fpB.A0 <- plogis(XfpA %*% fpApars.B0 + XfpA.offset)
    # False-positive prob of B, possibly dependent upon presence of A
    fpBpars <- fpBpars.A0 <- pars[(nOPA+nOPB+nG+nPA+nPB+nFPA+1):
                                  (nOPA+nOPB+nG+nPA+nPB+nFPA+nFPB)]
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



}

