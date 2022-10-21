#include <RcppArmadillo.h>
#include <float.h>
#include "distprob.h"
#include "pifun.h"
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gdistremoval(arma::vec beta, arma::uvec n_param, arma::vec yDistance,
    arma::vec yRemoval, arma::mat ysum, int mixture, std::string keyfun,
    arma::mat Xlam, arma::vec A, arma::mat Xphi, arma::mat Xrem,
    arma::mat Xdist, arma::vec db, arma::mat a, arma::mat u, arma::vec w,
    arma::uvec pl, int K, arma::uvec Kmin, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  int M = Xlam.n_rows;
  int T = Xphi.n_rows / M;
  int Rdist = yDistance.size() / M;
  unsigned Jdist = Rdist / T;
  int Rrem = yRemoval.size() / M;
  unsigned Jrem = Rrem / T;

  //Abundance
  const vec lambda = exp(Xlam * beta_sub(beta, n_param, 0)) % A;
  double log_alpha = beta_sub(beta, n_param, 1)(0); //length 1 vector

  //Availability
  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit(Xphi * beta_sub(beta, n_param, 2));
  }

  //Distance sampling detection
  vec dist_param(M*T);
  if(keyfun != "uniform"){
    dist_param = exp(Xdist * beta_sub(beta, n_param, 3));
  }
  double scale = exp(beta_sub(beta, n_param, 4)(0));

  //Removal sampling detection
  vec remP(M * T * Jrem);
  remP = inv_logit(Xrem * beta_sub(beta, n_param, 5));

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){
    double site_lp = 0.0;
    vec f = zeros(K+1);
    for (int k=Kmin(i); k<(K+1); k++){
      f(k) = N_density(mixture, k, lambda(i), log_alpha);
    }

    int t_ind = i * T;
    int yd_ind = i * T * Jdist;
    int yr_ind = i * T * Jrem;

    vec yd_sub(Jdist);
    vec yr_sub(Jrem);
    vec cpd(Jdist);
    vec cpr(Jrem);
    double pdist;
    double prem;
    vec g = ones(K+1);

    for(int t=0; t<T; t++){
      int yd_stop = yd_ind + Jdist - 1;
      yd_sub = yDistance.subvec(yd_ind, yd_stop);

      int yr_stop = yr_ind + Jrem - 1;
      yr_sub = yRemoval.subvec(yr_ind, yr_stop);

      uvec nmd = find_finite(yd_sub);
      uvec nmr = find_finite(yr_sub);

      if((nmd.size() == Jdist) && (nmr.size() == Jrem)){

        cpd = distprob(keyfun, dist_param(t_ind), scale, "point",
                           db, w, a.row(i)) % u.col(i);
        pdist = sum(cpd);

        cpr = removalPiFun(remP.subvec(yr_ind, yr_stop), pl);
        prem = sum(cpr);

        //if(pdist == 0 | prem == 0){
        //  site_lp += R_PosInf;
        //}

        site_lp += dmultinom(yd_sub, cpd/pdist);
        site_lp += dmultinom(yr_sub, cpr/prem);

        for (int k=Kmin(i); k<(K+1); k++){
          g(k) *= Rf_dbinom(ysum(i,t), k, pdist*prem*phi(t_ind), false);
        }

      }

      t_ind += 1;
      yd_ind += Jdist;
      yr_ind += Jrem;
    }

    site_lp += log(sum(f % g));

    loglik += site_lp + DBL_MIN;

  }

  return -loglik;

}
