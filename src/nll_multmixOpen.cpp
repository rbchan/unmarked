#include "nll_multmixOpen.h"

using namespace Rcpp;
using namespace arma;

SEXP nll_multmixOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_,
    SEXP Xp_, SEXP Xiota_,
    SEXP beta_, SEXP beta_ind_,
    SEXP Xlam_offset_,
    SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP Xiota_offset_,
    SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_,
    SEXP first_, SEXP last_, SEXP first1_,
    SEXP M_, SEXP T_, SEXP J_, SEXP R_, SEXP delta_, SEXP dynamics_,
    SEXP fix_, SEXP go_dims_, SEXP immigration_,
    SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP pi_fun_,
    SEXP lfac_k_, SEXP kmyt_, SEXP lfac_kmyt_, SEXP fin_) {

  //Indices
  int lk = as<int>(lk_);
  Rcpp::IntegerVector N = seq_len(lk)-1;
  int M = as<int>(M_);
  int T = as<int>(T_);
  int J = as<int>(J_);
  int R = as<int>(R_);
  ucube y = as<ucube>(y_);
  imat yt = as<imat>(yt_);
  Rcpp::IntegerVector first(first_);
  Rcpp::IntegerVector last(last_);
  int first1 = as<int>(first1_);
  arma::imat ytna = as<arma::imat>(ytna_); // y[i,,t] are all NA
  ucube yna = as<ucube>(yna_);
  arma::imat delta = as<arma::imat>(delta_);

  vec lfac_k = as<vec>(lfac_k_);
  cube lfac_kmyt = as<cube>(lfac_kmyt_);
  cube kmyt = as<cube>(kmyt_);
  icube fin = as<icube>(fin_);
  imat I = as<arma::imat>(I_);
  imat I1 = as<arma::imat>(I1_);
  List Ib(Ib_);
  List Ip(Ip_);
  int nrI = I.n_rows;
  int nrI1 = I1.n_rows;

  //Covariate matrices
  mat Xlam = as<mat>(Xlam_);
  mat Xgam = as<mat>(Xgam_);
  mat Xom = as<mat>(Xom_);
  mat Xp = as<mat>(Xp_);
  mat Xiota = as<mat>(Xiota_);

  //Offsets
  colvec Xlam_offset = as<colvec>(Xlam_offset_);
  colvec Xgam_offset = as<colvec>(Xgam_offset_);
  colvec Xom_offset = as<colvec>(Xom_offset_);
  colvec Xp_offset = as<colvec>(Xp_offset_);
  colvec Xiota_offset = as<colvec>(Xiota_offset_);

  //Model types
  std::string mixture = as<std::string>(mixture_);
  std::string pi_fun = as<std::string>(pi_fun_);
  std::string dynamics = as<std::string>(dynamics_);
  std::string fix = as<std::string>(fix_);
  std::string go_dims = as<std::string>(go_dims_);
  bool immigration = as<bool>(immigration_);

  //Parameters
  vec beta = as<vec>(beta_);
  umat bi = as<umat>(beta_ind_);

  //Lambda
  vec beta_lam = beta.subvec(bi(0,0), bi(0,1));
  vec lam = exp(Xlam*beta_lam + Xlam_offset);

  //Get 2nd abundance dist parameter set up if necessary
  double alpha = 0.0;
  double psi = 0.0;
  if(mixture=="NB"){
    alpha = exp(beta(bi(5,0)));
  } else if(mixture=="ZIP"){
    psi = 1.0 / (1.0 + exp(-1 * beta(bi(5,0))));
  }

  //Omega
  mat omv = ones<colvec>(M*(T-1));
  if((fix != "omega") && (dynamics != "trend")) {
    vec beta_om = beta.subvec(bi(2,0), bi(2,1));
    if((dynamics == "ricker")  || (dynamics == "gompertz")){
        omv = exp(Xom*beta_om + Xom_offset);
    } else if((dynamics == "constant")  || (dynamics == "autoreg") ||
              (dynamics == "notrend")){
        omv = 1.0/(1.0+exp(-1*(Xom*beta_om + Xom_offset)));
    }
  }
  omv.reshape(T-1, M);
  mat om = trans(omv);

  //Gamma
  mat gam = zeros(M, T-1);
  if(dynamics == "notrend") {
    mat lamMat = repmat(lam, 1, T-1);
    gam = (1-om) % lamMat;
  } else if (fix != "gamma") {
    vec beta_gam = beta.subvec(bi(1,0), bi(1,1));
    mat gamv = exp(Xgam*beta_gam + Xgam_offset);
    gamv.reshape(T-1, M);
    gam = trans(gamv);
  }

  //Detection prob
  vec beta_p = beta.subvec(bi(3,0), bi(3,1));
  vec pv = inv_logit(Xp * beta_p + Xp_offset);
  //Transform into cube (J x M x T)
  cube p((const double*) pv.begin(), R, M, T);

  //Immigration
  mat iota = zeros(M,T-1);
  if(immigration){
    vec beta_iota = beta.subvec(bi(4,0), bi(4,1));
    mat iotav = exp(Xiota*beta_iota + Xiota_offset);
    iotav.reshape(T-1, M);
    iota = trans(iotav);
  }

  // initialize
  double ll=0.0;
  double ll_i=0.0;
  int first_i=0;
  int last_i=0;
  colvec g_star = ones(lk);
  mat g3 = zeros(lk, lk);
  mat g3_d = zeros(lk, lk);
  colvec g1_t = zeros(lk);
  colvec g1_t_star = zeros(lk);
  colvec g1 = zeros(lk);
  colvec g1_star = zeros(lk);
  cube g3_t = zeros(lk,lk,T-1);
  mat ky_slice(lk, T);
  umat y_slice(J,T);

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
	      tp1(g3_t.slice(t), nrI, nrI1, N, I, I1, Ib, Ip,
            gam(first1,t), om(first1,t));
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
    ll_i=0.0;
    first_i = first[i];
    last_i = last[i];

    ky_slice = kmyt.slice(i); //using subcube causes segfaults
    y_slice = y.slice(i);

    //Calculate g_star
    g_star.ones();
    if(last_i > first_i) {
      // loop over time periods in reverse order, up to second occasion
      for(int t=last_i; t>first_i; t--) {
        g1_t.zeros();
        //Skip site i time t if all NA
        if(ytna(i,t)==1) {
	        continue;
	      }

        //Detection
        uvec na_sub = yna.slice(i).col(t);

        //vec ysub = y.subcube(span(i), span(t), span());
        uvec ysub = y_slice.col(t);
        vec lkmyt_sub = lfac_kmyt.subcube(span(i), span(t), span());
        vec psub = p.subcube(span(), span(i), span(t));
        vec cp = piFun(psub, pi_fun);

        ivec fin_sub = fin.subcube(span(i), span(t), span());
        vec kmyt_sub = ky_slice.col(t);

        if(any(na_sub)){
          uvec ids = find(na_sub!=1);
          ysub = ysub.elem(ids);
          cp = cp.elem(ids);
        }

        double ycp = sum(ysub % log(cp));
        double cpJ = 1 - sum(cp);

        g1_t = lfac_k - lkmyt_sub + ycp + log(cpJ) * kmyt_sub;
        g1_t = exp(g1_t);
        g1_t = g1_t % fin_sub; //set to 0 when K impossible

	      g1_t_star = g1_t % g_star;

        // computes transition probs for g3
        if(go_dims == "matrix") {
	        g3.zeros();
	        if(dynamics=="constant" || dynamics=="notrend") {
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
        }
        else if(go_dims == "rowvec") {
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
        }
        else {
	        g_star = g3 * g1_t_star;
        }
      }
    }

    //Calculate g1
    int delta_i0 = delta(i,0);
    g1.zeros();

    uvec na_sub = yna.slice(i).col(first_i);

    //vec ysub = y.subcube(span(i), span(first_i), span());
    uvec ysub = y_slice.col(first_i);
    vec lkmyt_sub = lfac_kmyt.subcube(span(i), span(first_i), span());
    vec psub = p.subcube(span(), span(i), span(first_i));
    vec cp = piFun(psub, pi_fun);

    ivec fin_sub = fin.subcube(span(i), span(first_i), span());
    vec kmyt_sub = ky_slice.col(first_i);

    if(any(na_sub)){
      uvec ids = find(na_sub!=1);
      ysub = ysub.elem(ids);
      cp = cp.elem(ids);
    }

    double ycp = sum(ysub % log(cp));
    double cpJ = 1 - sum(cp);

    g1 = lfac_k - lkmyt_sub + ycp + log(cpJ) * kmyt_sub;
    g1 = exp(g1);
    g1 = g1 % fin_sub; //set to 0 when K impossible

    if(delta_i0>1){
	    g1_star = g1_t % g_star;
    }

    //Calculate g2
    //Nest k loops inside if statements to save time (?)
    colvec g2 = zeros(lk);
    if(mixture=="P"){
      for(int k=0; k<lk; k++){
        g2(k) = Rf_dpois(k, lam(i), false);
      }
    }
    else if(mixture=="NB") {
      for(int k=0; k<lk; k++){
        g2(k) = dnbinom_mu(k, alpha, lam(i), false);
      }
    }
    else if(mixture=="ZIP") {
      for(int k=0; k<lk; k++){
        g2(k) = dzip(k, lam(i), psi);
      }
    }

    //Combine g1, g2, gstar
    if(delta_i0==1) {
      vec g12star = g1 % g2 % g_star;
      ll_i += sum(g12star);
    }
    else if(delta_i0>1) {
      g3_d = g3;
      for(int d=0; d<delta_i0; d++) {
        g3_d = g3_d * g3;
      }
      g_star = g3_d * g1_star;
      vec g2star = g2 % g_star;
      ll_i += sum(g2star);
    }

    ll += log(ll_i + DBL_MIN);

  }

  return wrap(-ll);
}
