#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec nll_occuMulti_loglik(Rcpp::IntegerVector fStart, Rcpp::IntegerVector fStop,
    arma::sp_mat dmF, Rcpp::List dmOcc,
    arma::colvec beta, Rcpp::List dmDet,
    Rcpp::IntegerVector dStart, Rcpp::IntegerVector dStop,
    arma::mat y, Rcpp::IntegerVector yStart, Rcpp::IntegerVector yStop,
    arma::mat Iy0, arma::mat z, Rcpp::LogicalVector fixed0){

  int nF = dmF.n_rows; //dmF is already transposed
  int S = y.n_cols;
  int J = y.n_rows;
  int N = yStart.size();

  //psi calculation
  int index = 0;
  mat f(N, nF);
  for(int i = 0; i < nF; i++){
    if(fixed0(i)){
      f.col(i) = zeros(N);
    } else {
      mat X = as<mat>(dmOcc[index]);
      f.col(i) = X * beta.subvec(fStart[index], fStop[index]);
      index += 1;
    }
  }

  mat psi = exp( f * dmF );
  for(unsigned int i = 0; i < psi.n_rows; i++){
    psi.row(i) = psi.row(i) / sum( psi.row(i) );
  }

  //p calculation
  mat p(J, S);
  for(int i = 0; i < S; i++){
    mat X = as<mat>(dmDet[i]);
    colvec beta_sub = beta.subvec(dStart[i], dStop[i]);
    p.col(i) = 1.0 / ( 1.0 + exp( -1.0 * (X * beta_sub) ));
  }

  //Probability of detection history
  int M = z.n_rows;
  vec logLik(N);

  for(int i = 0; i < N; i++){

    mat ysub = y.rows(yStart[i], yStop[i]);
    mat psub = p.rows(yStart[i], yStop[i]);
    rowvec cdp(S);
    for(int j = 0; j < S; j++){
     cdp(j) = exp( sum( ysub.col(j) % log(psub.col(j)) ) +
                   sum( (1 - ysub.col(j)) % log( 1 - psub.col(j)) ) );
    }

    rowvec prdProbY(M);
    for(int j = 0; j < M; j++){
      prdProbY(j) = prod( z.row(j) % cdp + (1 - z.row(j)) % Iy0.row(i) );
    }

    logLik(i) = log( sum( psi.row(i) % prdProbY ) );

  }

  return logLik;
}

// [[Rcpp::export]]
double nll_occuMulti(Rcpp::IntegerVector fStart, Rcpp::IntegerVector fStop,
    arma::sp_mat dmF, Rcpp::List dmOcc,
    arma::colvec beta, Rcpp::List dmDet,
    Rcpp::IntegerVector dStart, Rcpp::IntegerVector dStop,
    arma::mat y, Rcpp::IntegerVector yStart, Rcpp::IntegerVector yStop,
    arma::mat Iy0, arma::mat z, Rcpp::LogicalVector fixed0, double penalty){

  vec logLik = nll_occuMulti_loglik(fStart, fStop, dmF, dmOcc, beta, dmDet,
      dStart, dStop, y, yStart, yStop, Iy0, z, fixed0);

  double pen = penalty * 0.5 * accu(pow(beta, 2));
  return -1.0 * (sum(logLik) - pen);
}
