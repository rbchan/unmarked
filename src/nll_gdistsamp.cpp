#include <RcppArmadillo.h>
#include <float.h>
#include "distprob.h"
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gdistsamp(arma::vec beta, arma::uvec n_param, arma::vec y,
    int mixture, std::string keyfun, std::string survey,
    arma::mat Xlam, arma::vec Xlam_offset, arma::vec A, arma::mat Xphi,
    arma::vec Xphi_offset, arma::mat Xdet, arma::vec Xdet_offset, arma::vec db,
    arma::mat a, arma::mat u, arma::vec w, int K, arma::uvec Kmin, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  int M = Xlam.n_rows;
  int T = Xphi.n_rows / M;
  unsigned J = db.size() - 1;

  //Abundance
  const vec lambda = exp(Xlam * beta_sub(beta, n_param, 0) + Xlam_offset) % A;
  double log_alpha = beta_sub(beta, n_param, 4)(0); //length 1 vector

  //Availability
  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit(Xphi * beta_sub(beta, n_param, 1) + Xphi_offset);
  }

  //Detection
  vec det_param(M*T);
  if(keyfun != "uniform"){
    det_param = exp(Xdet * beta_sub(beta, n_param, 2) + Xdet_offset);
  }
  double scale = exp(beta_sub(beta, n_param, 3)(0));

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){

    vec f = zeros(K+1); // should this be ones(K+1);
    for (int k=Kmin(i); k<(K+1); k++){
      f(k) = N_density(mixture, k, lambda(i), log_alpha);
    }

    int t_ind = i * T;
    int y_ind = i * T * J;

    vec y_sub(J);
    vec y_all(J+1);

    vec cp(J);
    vec cp_all(J+1);
    double ptotal;

    vec g = zeros(K+1);

    for(int t=0; t<T; t++){
      int y_stop = y_ind + J - 1;
      y_sub = y.subvec(y_ind, y_stop);
      y_all.subvec(0, (J-1)) = y_sub;

      uvec nm = find_finite(y_sub);

      if((nm.size() == J)){

        cp = distprob(keyfun, det_param(t_ind), scale, survey, db,
                      w, a.row(i));
        cp = cp % u.col(i) * phi(t_ind);
        ptotal = sum(cp);

        cp_all.subvec(0, (J-1)) = cp;
        cp_all(J) = 1 - ptotal;

        for (int k=Kmin(i); k<(K+1); k++){
          y_all(J) = k - sum(y_sub);
          g(k) += dmultinom(y_all, cp_all);
        }

      }

      t_ind += 1;
      y_ind += J;
    }

    loglik += log(sum(f % exp(g)));

  }

  return -loglik;

}
