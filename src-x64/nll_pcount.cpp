#include "nll_pcount.h"
#include "distr.h"

using namespace Rcpp ;

SEXP nll_pcount( SEXP yR, SEXP Xr, SEXP Vr, SEXP beta_lamR, SEXP beta_pR, SEXP log_alphaR, SEXP X_offsetR, SEXP V_offsetR, SEXP naMatR, SEXP lkR, SEXP mixtureR ) {
  arma::imat y = as<arma::imat>(yR);
  arma::mat X = as<arma::mat>(Xr);
  arma::mat V = as<arma::mat>(Vr);
  arma::colvec beta_lam = as<arma::colvec>(beta_lamR);
  arma::colvec beta_p = as<arma::colvec>(beta_pR);
  double log_alpha = as<double>(log_alphaR);
  arma::colvec X_offset = as<arma::colvec>(X_offsetR);
  arma::colvec V_offset = as<arma::colvec>(V_offsetR);
  Rcpp::LogicalMatrix naMat(naMatR);
  std::string mixture = as<std::string>(mixtureR);
  double alpha = 0.0;
  if(mixture=="NB")
    alpha = exp(log_alpha);
  else if(mixture=="ZIP")
    alpha = 1.0/(1.0+exp(-log_alpha));
  int lk = as<int>(lkR);
  int R = y.n_rows;
  int J = y.n_cols;
  arma::colvec lam = exp(X*beta_lam + X_offset);
  arma::colvec logit_p = V*beta_p + V_offset;
  arma::mat pv = 1.0/(1.0+exp(-logit_p));
  pv.reshape(J,R);
  arma::mat p = trans(pv);
  double L_i=0.0;
  double ll=0.0;
  arma::colvec f(lk);
  f.zeros();
  arma::colvec g(lk);
  for(int i=0; i<R; i++) {
    g.zeros();
    L_i = 0.0;
    for(int k=0; k<lk; k++) {
      if(mixture=="P")
	f(k) = Rf_dpois(k, lam(i), false);
      else if(mixture=="NB")
        f(k) = dnbinom_mu(k, alpha, lam(i), false);
      else if(mixture=="ZIP")
	f(k) = dzip(k, lam(i), alpha);
      for(int j=0; j<J; j++) {
	// if(k >= y(i,j))
	if(!naMat(i,j))
	  g(k) += Rf_dbinom(y(i,j), k, p(i,j), true);
	// else
	// g(k) = 0.0;
      }
      g(k) = exp(g(k));
      L_i += f(k) * g(k);
    }
    ll += log(L_i + DOUBLE_XMIN);
  }
  return wrap(-ll);

}
