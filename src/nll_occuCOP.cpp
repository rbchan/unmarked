#include <RcppArmadillo.h>
#include <float.h>

using namespace Rcpp ;

// [[Rcpp::export]]
double nll_occuCOP(arma::icolvec y, arma::icolvec L,
    arma::mat Xpsi, arma::mat Xlambda,
    arma::colvec beta_psi, arma::colvec beta_lambda,
    Rcpp::LogicalVector removed) {

  // Number of sites M and obs J
  int M = Xpsi.n_rows;
  int J = y.n_elem / M;

  //Calculate psi back-transformed from logit
  arma::colvec psi = 1.0/(1.0+exp(-Xpsi*beta_psi));

  //Calculate lambda back-transformed from log
  arma::colvec lambda = exp(Xlambda*beta_lambda);

  double ll=0.0;
  int k=0; // counter
  // for each site i in 1:M
  for(int i=0; i<M; i++) {
    double iLambdaL=0.0; // init sum(lambda_ij * L_ij)
    double iN=0.0; // init sum(y) = total count of detec at site i
    int NbRemoved=0; // init count of the removed observations at site i
    for(int j=0; j<J; j++) {
      if(!removed(k)) {
        // If the observation is not removed from the analysis
        // (removed if there is a NA in y, L or in the relevant covariates for this site and obs)
        iLambdaL += lambda(k)*L(k);
        iN += y(k);
        NbRemoved += 1;
      }
      k++;
    }
    if ((!NbRemoved) < J) {
      if(iN>0) {
        ll += log(psi(i) * pow(iLambdaL, iN) / tgamma(iN + 1) * exp(-iLambdaL));
      } else {
        ll += log(psi(i) * exp(-iLambdaL) + 1-psi(i));
      }
    }
  }
  return -ll;
}
