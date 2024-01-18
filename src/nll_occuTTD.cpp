#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//namespace {

  vec get_Py(vec e_lamt, vec delta, vec naflag){

    //Remove NAs
    if(any(naflag)){
      uvec ids = find(naflag != 1);
      e_lamt = e_lamt.elem(ids);
      delta = delta.elem(ids);
    }
    int sum_delt = any(delta);

    vec out(2);
    out(0) = 1 - sum_delt;
    out(1) = prod(e_lamt);
    return(out);
  }

  mat get_phi(rowvec phi_raw){
    mat out(2,2);
    out.row(0) = phi_raw.subvec(0,1);
    out.row(1) = phi_raw.subvec(2,3);
    return(out);
  }

//}

// [[Rcpp::export]]
double nll_occuTTD( arma::vec beta, arma::vec y, arma::vec delta,
    arma::mat W, arma::mat V, arma::mat Xgam, arma::mat Xeps,
    arma::vec pind, arma::vec dind, arma::vec cind, arma::vec eind,
    std::string lpsi, std::string tdist,
    int N, int T, int J,
    arma::vec naflag){

  int ys = y.size();

  //Get psi values
  colvec raw_psi = W * beta.subvec(pind(0), pind(1));
  if(lpsi == "cloglog"){
    raw_psi = 1 - exp(-exp(raw_psi));
  } else {
    raw_psi = 1 / (1 + exp(-raw_psi));
  }
  const mat psi = join_rows(1-raw_psi, raw_psi);

  //Get lambda values
  const vec lam = exp(V * beta.subvec(dind(0), dind(1)));

  vec e_lamt(ys);
  if(tdist == "weibull"){
    double k = exp(beta(beta.size() - 1));
    for(int i=0; i<ys; i++){
      e_lamt(i) = pow(k*lam(i)*pow(lam(i)*y(i),(k-1)),delta(i)) *
          exp(-1*pow(lam(i)*y(i),k));
    }
  } else {
    for(int i=0; i<ys; i++){
      //exponential
      e_lamt(i) = pow(lam(i),delta(i)) * exp(-lam(i)*y(i));
    }
  }

  mat phi_raw(N*(T-1), 4);
  if(T > 1){
    colvec col = Xgam * beta.subvec(cind(0), cind(1));
    colvec ext = Xeps * beta.subvec(eind(0), eind(1));
    ext = 1 / (1 + exp(-ext));
    col = 1 / (1 + exp(-col));
    phi_raw.col(0) = 1-col;
    phi_raw.col(1) = col;
    phi_raw.col(2) = ext;
    phi_raw.col(3) = 1 - ext;
  }

  vec lik(N);
  int ystart = 0;
  int yend;
  int phi_index = 0;
  for (int n=0; n<N; n++){

    mat phi_prod = eye(2,2);
    if(T > 1){
      for(int t=0; t<(T-1); t++){
        yend = ystart + J - 1;
        vec Py_t = get_Py(e_lamt.subvec(ystart,yend),
            delta.subvec(ystart,yend),
            naflag.subvec(ystart,yend));
        mat D_ph = diagmat(Py_t);
        mat phi_t = get_phi(phi_raw.row(phi_index));
        phi_prod = phi_prod * (D_ph * phi_t);
        ystart += J;
        phi_index += 1;
      }
    }

    yend = ystart + J - 1;

    vec Py_T = get_Py(e_lamt.subvec(ystart,yend),
            delta.subvec(ystart,yend),
            naflag.subvec(ystart,yend));
    ystart += J;

    rowvec psi_phi = psi.row(n) * phi_prod;
    lik(n) = dot(psi_phi, Py_T);

  }

  return -sum(log(lik));

}
