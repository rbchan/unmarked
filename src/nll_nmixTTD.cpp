#include <RcppArmadillo.h>
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_nmixTTD(const arma::vec beta, const arma::vec y, const arma::vec delta,
    const arma::mat W, const arma::mat V, const arma::umat pinds,
    const std::string mixture, const std::string tdist, int N, int J, int K,
    const arma::vec naflag, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  //Get abundance lambda values
  const vec lamN = exp(W * beta.subvec(pinds(0,0), pinds(0,1)));

  //Get detection lambda values
  const vec lamP = exp(V * beta.subvec(pinds(1,0), pinds(1,1)));

  //Get alpha if negative binomial or zip
  double alpha = 1.0;
  if(mixture == "NB"){
    alpha = exp(beta(pinds(2,0)));
  } else if(mixture == "ZIP"){
    alpha = inv_logit(beta(pinds(2,0)));
  }

  //Get shape if weibull
  double shp = 1.0;
  if(tdist == "weibull"){
    shp = exp(beta(pinds(3,0)));
  }

  double loglik = 0.0;
  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int n=0; n<N; n++){

    vec pK(K+1);
    if(mixture == "P"){
      for (int k=0; k<(K+1); k++){
        pK(k) = R::dpois(k, lamN(n), 0);
      }
    } else if(mixture == "NB"){
      for (int k=0; k<(K+1); k++){
        pK(k) = dnbinom_mu(k, alpha, lamN(n), false);
      }
    } else if(mixture == "ZIP"){
      for (int k=0; k<(K+1); k++){
        pK(k) = dzip(k, lamN(n), alpha);
      }
    }

    int ystart = n * J;
    int yend = ystart + J - 1;
    vec nasub = naflag.subvec(ystart, yend);
    vec dsub = delta.subvec(ystart, yend);
    vec lamsub = lamP.subvec(ystart, yend);
    vec ysub = y.subvec(ystart, yend);

    if(any(nasub)){
      uvec ids = find(nasub != 1);
      dsub = dsub.elem(ids);
      lamsub = lamsub.elem(ids);
      ysub = ysub.elem(ids);
    }
    int ys = ysub.size();

    vec pY(K+1);
    pY(0) = 1 - any(dsub);

    if(tdist == "weibull"){
      for (int k=1; k<(K+1); k++){
        vec e_lamt(lamsub.size());
        for (int i=0; i<ys; i++){
          double lam = lamsub(i) * k;
          e_lamt(i) = pow(shp*lam*pow(lam*ysub(i),(shp-1)),dsub(i)) *
            exp(-1*pow(lam*ysub(i),shp));
        }
        pY(k) = prod(e_lamt);
      }
    } else {
      for (int k=1; k<(K+1); k++){
        vec e_lamt(lamsub.size());
        for(int i=0; i<ys; i++){
          //exponential
          double lam = lamsub(i) * k;
          e_lamt(i) = pow(lam,dsub(i)) * exp(-lam*ysub(i));
        }
        pY(k) = prod(e_lamt);
      }
    }

    loglik += log(dot(pK, pY));

  }

  return -loglik;

}
