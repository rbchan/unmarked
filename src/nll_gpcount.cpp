#include "nll_gpcount.h"
#include "distr.h"

using namespace Rcpp ;

SEXP nll_gpcount( SEXP y_, SEXP Xlam_, SEXP Xphi_, SEXP Xp_, SEXP beta_lam_, SEXP beta_phi_, SEXP beta_p_, SEXP log_alpha_, SEXP Xlam_offset_, SEXP Xphi_offset_, SEXP Xp_offset_, SEXP M_, SEXP mixture_, SEXP numPrimary_ ) {
  arma::mat ym = as<arma::mat>(y_); // Can't test for NAs if y is imat using is_finite??
  arma::mat yimax = max(ym, 1);
  arma::mat Xlam = as<arma::mat>(Xlam_);
  arma::mat Xphi = as<arma::mat>(Xphi_);
  arma::mat Xp = as<arma::mat>(Xp_); // needs to be a cube
  arma::colvec beta_lam = as<arma::colvec>(beta_lam_);
  arma::colvec beta_phi = as<arma::colvec>(beta_phi_);
  arma::colvec beta_p = as<arma::colvec>(beta_p_);
  double log_alpha = as<double>(log_alpha_);
  arma::colvec Xlam_offset = as<arma::colvec>(Xlam_offset_);
  arma::colvec Xphi_offset = as<arma::colvec>(Xphi_offset_);
  arma::colvec Xp_offset = as<arma::colvec>(Xp_offset_);
  int M = as<int>(M_);
  int lM = M+1;
  std::string mixture = as<std::string>(mixture_);
  int T = as<int>(numPrimary_);
  double alpha = 0.0;
  if(mixture=="NB")
    alpha = exp(log_alpha);
  else if(mixture=="ZIP")
    alpha = 1.0/(1.0+exp(-log_alpha));
  int R = ym.n_rows;
  int J = ym.n_cols / T;
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::colvec logit_phi = Xphi*beta_phi + Xphi_offset;
  arma::mat phiv = 1.0/(1.0+exp(-logit_phi));
  arma::colvec logit_p = Xp*beta_p + Xp_offset;
  arma::mat pv = 1.0/(1.0+exp(-logit_p));
  phiv.reshape(T,R);
  arma::mat phi = trans(phiv);
  pv.reshape(J*T,R);
  arma::mat pm = trans(pv);
  arma::cube y(R,J,T);
  arma::cube p(R,J,T);
  for(int q=0; q<(R*J*T); q++) {
    y(q) = ym(q);
    p(q) = pm(q);
  }
  //  Rprintf("made it %f \\n", 1.);
  double L=0.0, g=0.0, h=0.0;
  arma::vec f = arma::zeros<arma::vec>(lM);
  arma::mat gh = arma::zeros<arma::mat>(lM, lM);
  arma::vec ghi = arma::zeros<arma::vec>(lM);
  for(int i=0; i<R; i++) {
    //    Rprintf("log-like %f \\n", L);
    for(int m=0; m<lM; m++) {
      if(m < yimax(i))
	f(m) = log(0.0);
      else if(mixture=="P")
	f(m) = Rf_dpois(m, lam(i), true);
      else if(mixture=="NB")
        f(m) = dnbinom_mu(m, alpha, lam(i), true);
      else if(mixture=="ZIP")
	f(m) = log(dzip(m, lam(i), alpha));
    }
    ghi.zeros();
    for(int t=0; t<T; t++) {
      gh.zeros();
      for(int m=0; m<lM; m++) { // FIXME: increment from max(y[i,])
	for(int n=0; n<lM; n++) {
	  if((n > m) || (m < yimax(i))) {
	    gh(n,m) = log(0.0);
	    continue;
	  }
	  g = 0.0;
	  if(arma::is_finite(phi(i,t)))
	    g = Rf_dbinom(n, m, phi(i,t), true);
	  h = 0.0;
	  for(int j=0; j<J; j++) {
	    if(arma::is_finite(y(i,j,t))) { // true if not NA, NaN, +/-Inf
	      h += Rf_dbinom(y(i,j,t), n, p(i,j,t), true);
	    }
	  }
	  gh(n,m) = g + h;
	}
	ghi(m) += log(arma::accu(exp(gh.col(m)))); // sum over N(t)
      }
    }
    L -= log(arma::accu(exp(f + ghi))); // sum over M
  }
  return wrap(L);
}
