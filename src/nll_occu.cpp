#include "nll_occu.h"

using namespace Rcpp ;

SEXP nll_occu( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_psiR, SEXP beta_pR, SEXP ndR, SEXP knownOccR, SEXP navecR, SEXP X_offsetR, SEXP V_offsetR ) {
  arma::icolvec y = as<arma::icolvec>(yR);
  arma::mat X = as<arma::mat>(Xr);
  arma::mat V = as<arma::mat>(Vr);
  arma::colvec beta_psi = as<arma::colvec>(beta_psiR);
  arma::colvec beta_p = as<arma::colvec>(beta_pR);
  Rcpp::IntegerVector nd(ndR);
  Rcpp::LogicalVector knownOcc(knownOccR);
  Rcpp::LogicalVector navec(navecR);
  arma::colvec X_offset = as<arma::colvec>(X_offsetR);
  arma::colvec V_offset = as<arma::colvec>(V_offsetR);
  int R = X.n_rows;
  int J = y.n_elem / R;
  arma::colvec logit_psi = X*beta_psi + X_offset;
  arma::colvec logit_p = V*beta_p + V_offset;
  arma::colvec psi = 1.0/(1.0+exp(-logit_psi));
  arma::colvec p = 1.0/(1.0+exp(-logit_p));
  double ll=0.0;
  int k=0; // counter
  for(int i=0; i<R; i++) {
    double cp=1.0;
    for(int j=0; j<J; j++) {
      if(!navec(k))
	cp *= pow(p(k),y(k)) * pow(1-p(k), 1-y(k));
      k++;
    }
    if(knownOcc(i))
      psi(i) = 1.0;
    if(nd(i)==0)
      ll += log(cp * psi(i) + DOUBLE_XMIN);
    else if(nd(i)==1)
      ll += log(cp * psi(i) + (1-psi(i)) + DOUBLE_XMIN);
  }
  return wrap(-ll);
}

