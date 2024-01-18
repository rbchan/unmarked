#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//Build parameter from design matrix and beta
mat get_param(const List& dm_list, const vec& beta, const mat& ind){

  const int m = dm_list.size();
  const mat sz = as<mat>(dm_list[0]); //not ideal
  const int l = sz.n_rows;
  mat out(l,m);

  for(int i=0; i < m; i++){
    mat X = as<mat>(dm_list[i]);
    out.col(i) = X * beta.subvec(ind(i,0), ind(i,1));
  }

  return(out);
}

rowvec multinom_logit(const rowvec& lp){
  int sz = lp.size() + 1;
  rowvec exp_lp = exp(lp);
  double row_sum = sum(exp_lp) + 1;
  rowvec out(sz);
  out(0) = 1 / row_sum;
  for(int i=1; i<sz; i++){
    out(i) = exp_lp(i-1) / row_sum;
  }
  return(out);
}

//calculate final psi depending on parameterization
mat get_psi(const mat& lp, const std::string& prm){
  int R = lp.n_rows;

  if(prm == "multinomial"){
    //[ 1 - psi_1:psi_m, psi_1:psi_m ]
    mat out(R, lp.n_cols + 1);
    for (int i=0; i<R; i++){
     out.row(i) = multinom_logit(lp.row(i));
    }
    return(out);

  } else if(prm == "condbinom"){
    //[ 1-psi, psi * (1-R), psi * R ]
    mat psi_raw = 1 / ( 1 + exp(-lp));
    mat out(R,3);
    for(int i=0; i<R; i++){
      out(i,0) = 1 - psi_raw(i,0);
      out(i,1) = psi_raw(i,0) * (1 - psi_raw(i,1));
      out(i,2) = psi_raw(i,0) * psi_raw(i,1);
    }
    return(out);
  } else {
    stop("Invalid parameterization for get_psi");
  }
}

//calculate phi matrix for a site-year depending on parameterization
mat get_phi(int S, const rowvec& lp, const std::string& prm){

  mat out(S,S);

  if(prm == "multinomial"){
    int index = 0;
    for (int i=0; i<S; i++){ //row
      rowvec lp_row(S);
      for (int j=0; j<S; j++){ //col
        if(i == j){
          lp_row(j) = 1;
        } else {
          lp_row(j) = exp(lp(index));
          index += 1;
        }
      }
      out.row(i) = lp_row / sum(lp_row);
    }
    return(out);

  } else if(prm == "condbinom"){
    rowvec lp_logit = 1 / ( 1 + exp(-lp));
    for(int i=0; i<S; i++){
      out(i,0) = 1 - lp_logit(i);
      out(i,1) = lp_logit(i) * (1 - lp_logit(i+3));
      out(i,2) = lp_logit(i) * lp_logit(i+3);
    }
    return(out);
  } else {
    stop("Invalid parameterization passed to get_phi");
  }
}

mat get_sdp(int S, const rowvec& lp, const mat& guide,
    const std::string& prm){

  mat out = zeros(S,S);

  if(prm == "multinomial"){
    for (unsigned int i=0; i<lp.size(); i++){
      out( guide(i,0), guide(i,1) ) = exp(lp(i));
    }
    for (int s=0; s<S; s++){
      out(s,0) = 1;
      double row_sum = sum(out.row(s));
      for (int j=0; j<S; j++){
        out(s,j) = out(s,j) / row_sum;
      }
    }
    return(out);

  } else if(prm == "condbinom"){
    //input probs order is p_1, p_2, delta
    rowvec probs = 1 / ( 1 + exp(-lp));
    out(0,0) = 1;
    out(1,0) = 1 - probs(0);
    out(1,1) = probs(0);
    out(2,0) = 1 - probs(1);
    out(2,1) = probs(1) * (1 - probs(2));
    out(2,2) = probs(1) * probs(2);
    return(out);
  } else{
    stop("Invalid parameterization passed to get_sdp");
  }
}

vec get_ph(const int S, const rowvec& y, const mat& probs,
    const rowvec& navec, const mat& guide, const std::string prm){

  int J = probs.n_rows;

  vec out = ones(S);
  for(int j=0; j<J; j++){
    if(! navec(j)){
      rowvec prsub = probs.row(j);
      mat sdp = get_sdp(S, prsub, guide, prm);
      for(int s=0; s<S; s++){
        out(s) = out(s) * sdp(s, y(j));
      }
    }
  }
  return(out);
}

// [[Rcpp::export]]
double nll_occuMS(arma::vec beta, arma::mat y,
    Rcpp::List dm_state, Rcpp::List dm_phi, Rcpp::List dm_det,
    arma::mat sind, arma::mat pind, arma::mat dind, std::string prm,
    int S, int T, int J, int N,
    arma::mat naflag, arma::mat guide){

  //Get psi values
  const mat raw_psi = get_param(dm_state, beta, sind);
  const mat psi = get_psi(raw_psi, prm);

  //Get phi values
  mat raw_phi;
  if(T>1){
    raw_phi = get_param(dm_phi, beta, pind);
  }

  //Get p values
  const mat p = get_param(dm_det, beta, dind);

  vec lik = zeros(N);
  int pstart = 0;
  int phi_index = 0;
  int pend;
  int ystart;
  int yend;

  for(int n=0; n<N; n++){
    rowvec ysub = y.row(n);
    rowvec nasub = naflag.row(n);
    ystart = 0;
    mat phi_prod = eye(S,S);

    if(T>1){
      for(int t=0; t<(T-1); t++){
        pend = pstart + J - 1;
        yend = ystart + J - 1;

        vec ph_t = get_ph(S, ysub.subvec(ystart,yend),
            p.rows(span(pstart, pend)),
            nasub.subvec(ystart,yend), guide, prm);
        mat D_ph = diagmat(ph_t);
        mat phi_t = get_phi(S, raw_phi.row(phi_index), prm);
        phi_prod = phi_prod * (D_ph * phi_t);
        pstart += J;
        ystart += J;
        phi_index += 1;
      }
    }

    pend = pstart + J - 1;
    yend = ystart + J - 1;

    vec ph_T = get_ph(S, ysub.subvec(ystart,yend),
        p.rows(span(pstart, pend)),
        nasub.subvec(ystart,yend), guide, prm);
    pstart += J;

    rowvec psi_phi = psi.row(n) * phi_prod;
    lik(n) = dot(psi_phi, ph_T);

  }

  return -sum(log(lik));

}
