#include <RcppArmadillo.h>
#include <float.h>

using namespace Rcpp ;

// [[Rcpp::export]]
double nll_occu(arma::icolvec y, arma::mat X, arma::mat V,
    arma::colvec beta_psi, arma::colvec beta_p,
    Rcpp::IntegerVector nd, Rcpp::LogicalVector knownOcc, Rcpp::LogicalVector navec,
    arma::colvec X_offset, arma::colvec V_offset, std::string link_psi) {

  int R = X.n_rows;
  int J = y.n_elem / R;

  //Calculate psi
  arma::colvec psi_lp = X*beta_psi + X_offset;
  arma::colvec psi(R);
  if(link_psi == "cloglog"){
    psi = 1 - exp(-exp(psi_lp));
  } else {
    psi = 1.0/(1.0+exp(-psi_lp));
  }

  //Calculate p
  arma::colvec logit_p = V*beta_p + V_offset;
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
      ll += log(cp * psi(i) + DBL_MIN);
    else if(nd(i)==1)
      ll += log(cp * psi(i) + (1-psi(i)) + DBL_MIN);
  }
  return -ll;
}
