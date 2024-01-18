#include <RcppArmadillo.h>
#include "pifun.h"

using namespace Rcpp;
using namespace arma;

mat inv_logit_( mat inp ){
  return(1 / (1 + exp(-1 * inp)));
}

// [[Rcpp::export]]
double nll_multinomPois(arma::vec beta, std::string pi_fun,
    arma::mat Xlam, arma::vec Xlam_offset, arma::mat Xdet, arma::vec Xdet_offset,
    arma::vec y, arma::vec navec, int nP, int nAP){

  int M = Xlam.n_rows;
  vec lambda = exp( Xlam * beta.subvec(0, (nAP - 1) ) + Xlam_offset );

  int J = Xdet.n_rows / M;
  int R = y.size() / M;
  vec p = inv_logit_( Xdet * beta.subvec(nAP,(nP-1)) + Xdet_offset);

  int y_ind = 0;
  int p_ind = 0;

  mat ll = zeros(M,R);
  for (int m=0; m<M; m++){

    int y_stop = y_ind + R - 1;
    int p_stop = p_ind + J - 1;

    vec na_sub = navec.subvec(y_ind, y_stop);
    if( ! all(na_sub) ){

      vec pi_lam = piFun( p.subvec(p_ind, p_stop), pi_fun ) * lambda(m);

      for (int r=0; r<R; r++){
        if( ! na_sub(r) ){
           ll(m,r) = R::dpois(y(y_ind+r), pi_lam(r), 1);
        }
      }
    }

      y_ind += R;
      p_ind += J;
  }

  return -accu(ll);

}
