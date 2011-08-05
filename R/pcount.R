
#' Fit the N-mixture point count model

pcount <- function(formula, data, K, mixture = c("P", "NB"), starts,
	method = "BFGS", control = list(), se = TRUE)
{
	mixture <- match.arg(mixture)
	if(!is(data, "unmarkedFramePCount"))
		stop("Data is not an unmarkedFramePCount object.")

	designMats <- getDesign(data, formula)
	X <- designMats$X; V <- designMats$V; y <- designMats$y
        X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
        if (is.null(X.offset)) {
          X.offset <- rep(0, nrow(X))
        }
        if (is.null(V.offset)) {
          V.offset <- rep(0, nrow(V))
        }

	J <- ncol(y)
	M <- nrow(y)

	lamParms <- colnames(X)
	detParms <- colnames(V)
	nDP <- ncol(V)
	nAP <- ncol(X)

	if(missing(K)) K <- max(y, na.rm = TRUE) + 100
	if(K <= max(y, na.rm = TRUE))
		stop("specified K is too small. Try a value larger than any observation")
	k <- 0:K
	M <- nrow(y)
	J <- ncol(y)
	k.ik <- rep(k, M)
	k.ijk <- rep(k, M*J)

	nP <- nAP + nDP + ifelse(identical(mixture,"NB"),1,0)
	if(!missing(starts) && length(starts) != nP)
	   stop(paste("The number of starting values should be", nP))

	y.ij <- as.numeric(t(y))
	y.ijk <- rep(y.ij, each = K + 1)
	navec <- is.na(y.ijk)
	nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # I(no detection at site i)
	ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
	ijk.to.ikj <- with(ijk, order(i, k, j))

	nll <- function(parms){
		theta.i <- exp(X %*% parms[1 : nAP] + X.offset)
		p.ij <- plogis(V %*% parms[(nAP + 1) : (nAP + nDP)] + V.offset)
		theta.ik <- rep(theta.i, each = K + 1)
		p.ijk <- rep(p.ij, each = K + 1)

		bin.ijk <- dbinom(y.ijk,k.ijk,p.ijk)
		bin.ijk[which(is.na(bin.ijk))] <- 1
		bin.ik.mat <- matrix(bin.ijk[ijk.to.ikj], M * (K + 1), J, byrow = TRUE)
		g.ik <- rowProds(bin.ik.mat)

		if(identical(mixture,"P")) {
			f.ik <- dpois(k.ik,theta.ik)
		}
		else if (identical(mixture,"NB")){
			f.ik <- dnbinom(k.ik, mu = theta.ik, size = exp(parms[nP]))
		}
		dens.i.mat <- matrix(f.ik * g.ik, M, K + 1, byrow = TRUE)
		dens.i <- rowSums(dens.i.mat)  # sum over the K

		-sum(log(dens.i))
		}

	if(missing(starts)) starts <- rep(0, nP)
	fm <- optim(starts, nll, method=method, hessian = se, control = control)
	opt <- fm

	ests <- fm$par
	if(identical(mixture,"NB"))
		nbParm <- "alpha"
	else
		nbParm <- character(0)
	names(ests) <- c(lamParms, detParms, nbParm)
	if(se) {
		tryCatch(covMat <- solve(fm$hessian),
				error=function(x) stop(simpleError("Hessian is singular.  Try using fewer covariates.")))
	} else {
		covMat <- matrix(NA, nP, nP)
	}
	fmAIC <- 2 * fm$value + 2 * nP

	stateName <- "Abundance"

	stateEstimates <- unmarkedEstimate(name = stateName, short.name = "lam",
		estimates = ests[1:nAP], covMat = as.matrix(covMat[1:nAP,1:nAP]),
		invlink = "exp", invlinkGrad = "exp")

	detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
		estimates = ests[(nAP + 1) : (nAP + nDP)],
		covMat = as.matrix(covMat[(nAP + 1):(nAP + nDP), (nAP + 1):(nAP + nDP)]),
		invlink = "logistic", invlinkGrad = "logistic.grad")

	estimateList <- unmarkedEstimateList(list(state=stateEstimates,
		det=detEstimates))

	if(identical(mixture,"NB")) {
		estimateList@estimates$alpha <- unmarkedEstimate(name = "Dispersion",
			short.name = "alpha", estimates = ests[nP],
			covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
			invlinkGrad = "exp")
		}

	umfit <- new("unmarkedFitPCount", fitType = "pcount", call = match.call(),
		formula = formula, data = data, sitesRemoved = designMats$removed.sites,
		estimates = estimateList, AIC = fmAIC, opt = opt, negLogLike = fm$value,
		nllFun = nll, K = K, mixture = mixture)

	return(umfit)
}














pcountRcpp <- function(formula, data, K, mixture = c("P", "NB"), starts,
	method = "BFGS", control = list(), se = TRUE)
{
    mixture <- match.arg(mixture)
    if(!is(data, "unmarkedFramePCount"))
        stop("Data is not an unmarkedFramePCount object.")

    designMats <- unmarked:::getDesign(data, formula)
    X <- designMats$X; V <- designMats$V; y <- designMats$y
    X.offset <- designMats$X.offset; V.offset <- designMats$V.offset
    if (is.null(X.offset)) {
        X.offset <- rep(0, nrow(X))
    }
    if (is.null(V.offset)) {
        V.offset <- rep(0, nrow(V))
    }

    J <- ncol(y)
    M <- nrow(y)

    lamParms <- colnames(X)
    detParms <- colnames(V)
    nDP <- ncol(V)
    nAP <- ncol(X)

    if(missing(K)) K <- max(y, na.rm = TRUE) + 100
    if(K <= max(y, na.rm = TRUE))
        stop("specified K is too small. Try a value larger than any observation")
    k <- 0:K
    M <- nrow(y)
    J <- ncol(y)
    k.ik <- rep(k, M)
    k.ijk <- rep(k, M*J)

    nP <- nAP + nDP + ifelse(identical(mixture,"NB"),1,0)
    if(!missing(starts) && length(starts) != nP)
        stop(paste("The number of starting values should be", nP))

    y.ij <- as.numeric(t(y))
    y.ijk <- rep(y.ij, each = K + 1)
    navec <- is.na(y.ijk)
    nd <- ifelse(rowSums(y, na.rm=TRUE) == 0, 1, 0) # no detection
    ijk <- expand.grid(k = 0:K, j = 1:J, i = 1:M)
    ijk.to.ikj <- with(ijk, order(i, k, j))

    nll <- function(parms) {
        beta.lam <- parms[1:nAP]
        beta.p <- parms[(nAP+1):nP]
        nll.pcount(y, X, V, beta.lam, beta.p, K+1)
    }

    if(missing(starts)) starts <- rep(0, nP)
    fm <- optim(starts, nll, method=method, hessian = se, control = control)
    opt <- fm

    ests <- fm$par
    if(identical(mixture,"NB"))
        nbParm <- "alpha"
    else
        nbParm <- character(0)
    names(ests) <- c(lamParms, detParms, nbParm)
    if(se) {
        tryCatch(covMat <- solve(fm$hessian),
                 error=function(x) stop(simpleError("Hessian is singular.  Try using fewer covariates.")))
    } else {
        covMat <- matrix(NA, nP, nP)
    }
    fmAIC <- 2 * fm$value + 2 * nP

    stateName <- "Abundance"

    stateEstimates <- unmarkedEstimate(name = stateName,
        short.name = "lam", estimates = ests[1:nAP],
        covMat = as.matrix(covMat[1:nAP,1:nAP]),
        invlink = "exp", invlinkGrad = "exp")

    detEstimates <- unmarkedEstimate(name = "Detection", short.name = "p",
	estimates = ests[(nAP + 1) : (nAP + nDP)],
        covMat = as.matrix(covMat[(nAP + 1):(nAP + nDP),
                                  (nAP + 1):(nAP + nDP)]),
        invlink = "logistic", invlinkGrad = "logistic.grad")

    estimateList <- unmarked:::unmarkedEstimateList(list(
                                              state=stateEstimates,
                                              det=detEstimates))

    if(identical(mixture,"NB")) {
        estimateList@estimates$alpha <- unmarkedEstimate(
           name = "Dispersion",
           short.name = "alpha", estimates = ests[nP],
           covMat = as.matrix(covMat[nP, nP]), invlink = "exp",
           invlinkGrad = "exp")
    }

    umfit <- new("unmarkedFitPCount", fitType = "pcount",
                 call = match.call(), formula = formula, data = data,
                 sitesRemoved = designMats$removed.sites,
                 estimates = estimateList, AIC = fmAIC, opt = opt,
                 negLogLike = fm$value, nllFun = nll, K = K,
                 mixture = mixture)

    return(umfit)
}





















src <-  "
        arma::imat y = as<arma::imat>(yR);
        arma::mat X = as<arma::mat>(Xr);
        arma::mat V = as<arma::mat>(Vr);
        arma::rowvec beta_lam = as<arma::rowvec>(beta_lamR);
        arma::rowvec beta_p = as<arma::rowvec>(beta_pR);
        int lk = as<int>(lkR);
        int R = X.n_rows;
        int J = y.n_elem / R;
        arma::vec lam = exp(X*beta_lam);
        Rcpp::NumericMatrix pmat(R,J);
        arma::vec logit_p = V*beta_p;
        arma::vec p = 1.0/(1.0+exp(-logit_p));
        double ll=0.0;
        Rcpp::NumericVector f(lk);
        Rcpp::NumericVector g(lk);
        int z, zi;
        for(int i=0; i<R; i++) {
          zi=i*J;
          for(int k=0; k<lk; k++) {
            f(k) = Rf_dpois(k, lam(i), false);
            g(k) = 1.0;
            for(int j=0; j<J; j++) {
              z = zi + j;
//              if(k >= y(i,j))
                g(k) *= Rf_dbinom(y(i,j), k, p(z), false);
//              else
//                g(k) = 0.0;
              }
            }
          ll += log(sum(f*g));
          }
        return wrap(-ll);
        "


nll.pcount <- cxxfunction(signature(yR="numeric", Xr="matrix", Vr="matrix",
                          beta_lamR="numeric", beta_pR="numeric",
                          lkR="integer"),
                        body=src, plugin="RcppArmadillo")




  for(int k=1;k<=S;k++) {
    f(k) = pow(lambda, N(k)) * exp(-1.0*lambda - gammln(N(k)+1));
    }

  for(int i=1;i<=R;i++) {
      p(i) = 1.0/(1.0+exp(-1.0*(p0 + p1*x(i))));
      for(int k=1;k<=S;k++) {
          g(k) = 1.0;
          for(int j=1;j<=T;j++) {
              if(N(k)>=y(i,j))
                  g(k) *= exp(gammln(N(k)+1) - gammln(y(i,j)+1) -
                      gammln(N(k)-y(i,j)+1)) * pow(p(i), y(i,j)) *
                      pow(1.0-p(i), N(k)-y(i,j));
              else
                  g(k) = 0.0;
              }
          }
      nll -= log(sum(elem_prod(f, g)));
      }

z <- 0
for(i in 1:4) {
    z.i <- (i-1)*3
    for(j in 1:3) {
        z <- z.i+j
        cat(z, "\n")
    }
}
