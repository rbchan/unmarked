#include <RcppArmadillo.h>
#include <float.h>
#include "distr.h"
#include "tranprobs.h"


using namespace Rcpp ;


// [[Rcpp::export]]
double nll_pcountOpen(arma::imat ym, arma::mat Xlam, arma::mat Xgam, arma::mat Xom,
    arma::mat Xp, arma::mat Xiota, arma::colvec beta_lam, arma::colvec beta_gam,
    arma::colvec beta_om, arma::colvec beta_p, arma::colvec beta_iota,
    double log_alpha, arma::colvec Xlam_offset, arma::colvec Xgam_offset,
    arma::colvec Xom_offset, arma::colvec Xp_offset, arma::colvec Xiota_offset,
    arma::imat ytna, arma::imat ynam, int lk, std::string mixture,
    Rcpp::IntegerVector first, Rcpp::IntegerVector last, int M, int J, int T,
    arma::imat delta, std::string dynamics, std::string fix, std::string go_dims,
    bool immigration, arma::imat I, arma::imat I1, Rcpp::List Ib, Rcpp::List Ip) {

  Rcpp::IntegerVector N = seq_len(lk)-1;
  int nrI = I.n_rows;
  int nrI1 = I1.n_rows;
  double alpha=0.0, psi=0.0;
  if(mixture=="NB")
    alpha = exp(log_alpha);
  else if(mixture=="ZIP")
    psi = 1.0/(1.0+exp(-log_alpha));

  // linear predictors
  arma::colvec lam = exp(Xlam*beta_lam + Xlam_offset);
  arma::mat omv = arma::ones<arma::colvec>(M*(T-1));
  if((fix != "omega") && (dynamics != "trend")) {
    if((dynamics == "ricker")  || (dynamics == "gompertz"))
        omv = exp(Xom*beta_om + Xom_offset);
    else if((dynamics == "constant")  || (dynamics == "autoreg") || (dynamics == "notrend"))
        omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
  }
  omv.reshape(T-1, M);
  arma::mat om = arma::trans(omv);
  arma::mat gam = arma::zeros<arma::mat>(M, T-1);
  if(dynamics == "notrend") {
    arma::mat lamMat = arma::repmat(lam, 1, T-1);
    gam = (1-om) % lamMat;
  } else {
    if(fix != "gamma") {
      arma::mat gamv = exp(Xgam*beta_gam + Xgam_offset);
      gamv.reshape(T-1, M);
      gam = arma::trans(gamv);
    }
  }
  arma::mat pv = 1.0/(1.0+exp(-1*(Xp*beta_p + Xp_offset)));
  pv.reshape(J*T, M);
  arma::mat pm = trans(pv);
  //Immigration
  arma::mat iotav = arma::zeros<arma::colvec>(M*(T-1));
  if(immigration) {
    iotav = exp(Xiota*beta_iota + Xiota_offset);
  }
  iotav.reshape(T-1, M);
  arma::mat iota = arma::trans(iotav);
  // format matrices as cubes (shouldn't be done in likelihood)
  arma::icube y(M,J,T);
  arma::cube p(M,J,T);
  arma::icube yna(M,J,T);
  for(int q=0; q<(M*J*T); q++) {
    y(q) = ym(q);
    p(q) = pm(q);
    yna(q) = ynam(q);
  }
  // initialize
  double ll=0.0;
  double ll_i=0.0;
  int first_i=0;
  int last_i=0;
  int first1=0;
  arma::colvec g_star = arma::ones<arma::colvec>(lk);
  arma::mat g3 = arma::zeros<arma::mat>(lk, lk);
  arma::mat g3_d = arma::zeros<arma::mat>(lk, lk);
  arma::colvec g1_t = arma::zeros<arma::colvec>(lk);
  arma::colvec g1_t_star = arma::zeros<arma::colvec>(lk);
  arma::colvec g1 = arma::zeros<arma::colvec>(lk);
  arma::colvec g2 = arma::zeros<arma::colvec>(lk);
  arma::colvec g1_star = arma::zeros<arma::colvec>(lk);
  arma::cube g3_t = arma::zeros<arma::cube>(lk,lk,T-1);

  // shouldn't be done in likelihood
  for(int i=0; i<M; i++) {
    if(first[i]==1) {
      first1=i; // site with no missing values at t=1
      break;
    }
  }

  // compute g3 if there are no covariates of omega/gamma
  if(go_dims == "scalar") {
    if(dynamics=="constant" || dynamics=="notrend")
      tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,0), om(first1,0));
    else if(dynamics=="autoreg")
      tp2(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
    else if(dynamics=="trend")
      tp3(g3, lk, gam(first1,0), iota(first1,0));
    else if(dynamics=="ricker")
      tp4(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
    else if(dynamics=="gompertz")
      tp5(g3, lk, gam(first1,0), om(first1,0), iota(first1,0));
  } else if(go_dims == "rowvec") {
    for(int t=0; t<(T-1); t++) {
      if(ytna(first1,t)==1) { // FIXME: this is not generic!
	continue;
      }
      if(dynamics=="constant" || dynamics=="notrend") {
	tp1(g3_t.slice(t), nrI, nrI1, N, I, I1, Ib, Ip, gam(first1,t), om(first1,t));
      }
      else if(dynamics=="autoreg") {
	tp2(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
    }
      else if(dynamics=="trend")
	tp3(g3_t.slice(t), lk, gam(first1,t), iota(first1,t));
      else if(dynamics=="ricker")
	tp4(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
      else if(dynamics=="gompertz")
	tp5(g3_t.slice(t), lk, gam(first1,t), om(first1,t), iota(first1,t));
    }
  }

  // loop over sites
  for(int i=0; i<M; i++) {
    first_i = first[i]-1; // remember 0=1st location in C++
    last_i = last[i]-1;
    g_star.ones();
    if(last_i > first_i) {
      // loop over time periods in reverse order, up to second occasion
      for(int t=last_i; t>first_i; t--) {
	if(ytna(i,t)==1) {
	  continue; //
	}
	g1_t.zeros();
	// loop over possible value of N at time t
	for(int k=0; k<lk; k++) {
	  for(int j=0; j<J; j++) {
	    if(yna(i,j,t)==0) {
	      g1_t(k) += Rf_dbinom(y(i,j,t), k, p(i,j,t), true);
	    }
	  }
	  g1_t(k) = exp(g1_t(k));
	  g1_t_star(k) = g1_t(k) * g_star(k);
	}
	// computes transition probs for g3
	if(go_dims == "matrix") {
	  g3.zeros();
	  if(dynamics=="constant" || dynamics=="notrend") {
	    //	    tp1(g3, lk, gam(i,t-1), om(i,t-1));
	    tp1(g3, nrI, nrI1, N, I, I1, Ib, Ip, gam(i,t-1), om(i,t-1));
	  }
	  else if(dynamics=="autoreg") {
	    tp2(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	  }
	  else if(dynamics=="trend")
	    tp3(g3, lk, gam(i,t-1), iota(i,t-1));
	  else if(dynamics=="ricker")
	    tp4(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	  else if(dynamics=="gompertz")
	    tp5(g3, lk, gam(i,t-1), om(i,t-1), iota(i,t-1));
	} else if(go_dims == "rowvec") {
	  g3 = g3_t.slice(t-1);
	}
	int delta_it = delta(i,t);
	// matrix multiply transition probs over time gaps
	if(delta_it>1) {
	  g3_d = g3;
	  for(int d=1; d<delta_it; d++) {
	    g3_d = g3_d * g3;
	    }
	  g_star = g3_d * g1_t_star;
	} else
	  g_star = g3 * g1_t_star;
      }
    }
    ll_i=0.0;
    int delta_i0 = delta(i,0);
    g1.zeros();
    for(int k=0; k<lk; k++) { // loop over possible values of N
      for(int j=0; j<J; j++) {
	if(yna(i,j,first_i)==0) {
	  g1(k) += Rf_dbinom(y(i,j,first_i), k, p(i,j,first_i), true);
	}
      }
      g1(k) = exp(g1(k));
      if(delta_i0>1)
	g1_star(k) = g1(k) * g_star(k);
      if(mixture=="P")
	g2(k) = Rf_dpois(k, lam(i), false);
      else if(mixture=="NB")
        g2(k) = dnbinom_mu(k, alpha, lam(i), false);
      else if(mixture=="ZIP")
	g2(k) = dzip(k, lam(i), psi);
      if(delta_i0==1)
	ll_i += g1(k) * g2(k) * g_star(k);
    }
    if(delta_i0>1) {
      g3_d = g3;
      for(int d=0; d<delta_i0; d++) {
	g3_d = g3_d * g3;
      }
      g_star = g3_d * g1_star;
      ll_i = arma::dot(g2, g_star);
    }
    ll += log(ll_i + DBL_MIN);
  }
  return -ll;
}
