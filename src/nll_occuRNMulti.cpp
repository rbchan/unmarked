#include <RcppArmadillo.h>
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//Calculate part of abundance likelihood to integrate over (abun f[g] * detect [g])
double calc_fg(int k, double lam, int J, arma::vec r, arma::rowvec y,
               arma::irowvec miss){
  double f = R::dpois(k, lam, 0);
  double p;
  double g = 0;
  for (int j = 0; j<J; j++){
    if(miss(j)) continue;
    p = 1 - pow(1 - r(j), k);
    g += R::dbinom(y(j), 1, p, 1);
  }
  return f * exp(g);
}

//Calculate abundance likelihood for subordinate species
double lik_subord_abun(int Kmin, int Kmax, double lam, int J, arma::vec r,
    arma::rowvec y, arma::irowvec miss){

  double lik = 0.0;
  lam = exp(lam);
  for (int k = Kmin; k < (Kmax+1); k++){
    lik += calc_fg(k, lam, J, r, y, miss);
  }

  return lik;
}

//Calculate occupancy likelihood for subordinate species
double lik_subord_occ(int Kmin, double psi, int J, arma::vec p, arma::rowvec y,
                      arma::irowvec miss){
  double g = 1;
  psi = inv_logit(psi);
  for (int j = 0; j<J; j++){
    if(miss(j)) continue;
    g *= R::dbinom(y(j), 1, p(j), 0);
  }
  return g * psi + (1 - Kmin) * (1 - psi);
}

// [[Rcpp::export]]
double nll_occuRNMulti(arma::vec beta, arma::mat state_ind, arma::mat det_ind, int S,
    arma::ivec modOcc, Rcpp::List ylist, Rcpp::List Xlist, Rcpp::List Vlist,
    arma::imat dep, arma::ivec K, arma::imat Kmin,
    arma::imat miss, arma::ivec site_miss, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  double nll = 0;

  mat y1 = as<mat>(ylist[0]);
  int M = y1.n_rows;
  int J = y1.n_cols;

  //Parameters for dominant species
  mat X1 = as<mat>(Xlist[0]);
  vec lam1 = exp(X1 * beta.subvec(state_ind(0, 0), state_ind(0, 1)));

  mat V1 = as<mat>(Vlist[0]);
  vec r1 = inv_logit(V1 * beta.subvec(det_ind(0, 0), det_ind(0,1)));

  //Parameters for 2nd species
  mat y2 = as<mat>(ylist[1]);

  int sind_row = 1;

  vec lam2;
  vec lp2, lp2_sp1;
  if(dep(1,0)){
    List Xlist_sp2 = Xlist[1];
    mat X2 = as<mat>(Xlist_sp2[0]);
    lp2 = X2 * beta.subvec(state_ind(1,0), state_ind(1,1));
    mat X2_sp1 = as<mat>(Xlist_sp2[1]);
    lp2_sp1 = X2_sp1 * beta.subvec(state_ind(2,0), state_ind(2,1));
    sind_row += 2;
  } else {
    mat X2 = as<mat>(Xlist[1]);
    lam2 = exp(X2 * beta.subvec(state_ind(1,0), state_ind(1,1)));
    sind_row += 1;
  }

  mat V2 = as<mat>(Vlist[1]);
  vec r2 = inv_logit(V2 * beta.subvec(det_ind(1, 0), det_ind(1,1)));

  //Parameters for 3rd species (if needed)
  mat y3;
  vec r3;
  vec lam3, lp3, lp3_sp1, lp3_sp2;

  //Very janky...
  if(S > 2){
    y3 = as<mat>(ylist[2]);

    if(dep(2,0) | dep(2,1)){
      List Xlist_sp3 = Xlist[2];
      mat X3 = as<mat>(Xlist_sp3[0]);
      lp3 = X3 * beta.subvec(state_ind(sind_row, 0), state_ind(sind_row, 1));
      sind_row += 1;

      int mat_idx = 1;
      if(dep(2,0)){
        mat X3_sp1 = as<mat>(Xlist_sp3[mat_idx]);
        mat_idx += 1;
        lp3_sp1 = X3_sp1 * beta.subvec(state_ind(sind_row, 0), state_ind(sind_row, 1));
        sind_row += 1;
      }
      if(dep(2,1)){
        mat X3_sp2 = as<mat>(Xlist_sp3[mat_idx]);
        lp3_sp2 = X3_sp2 * beta.subvec(state_ind(sind_row, 0), state_ind(sind_row, 1));
      }
    } else {
      // Unused at the moment I think
      mat X3 = as<mat>(Xlist[2]);
      lam3 = exp(X3 * beta.subvec(state_ind(sind_row, 0), state_ind(sind_row, 1)));
    }

    mat V3 = as<mat>(Vlist[2]);
    r3 = inv_logit(V3 * beta.subvec(det_ind(2, 0), det_ind(2,1)));
  }

  #pragma omp parallel for reduction(-: nll) if(threads > 1)
  for (int m = 0; m<M; m++){
    double fg1, fg2, lik_sp2, lik_sp3;

    if(site_miss(m)) continue;

    //Indices for detection probability
    int r_st = m * J;
    int r_end = r_st + J - 1;
    vec r1_sub = r1.subvec(r_st, r_end);
    vec r2_sub = r2.subvec(r_st, r_end);

    vec r3_sub;
    if(S==3){
      r3_sub = r3.subvec(r_st, r_end);
    }

    double lik_m = 0.0;

    double par2_m, par3_m; // state parameters for species 2 and 3

    //This can probably be hugely simplified by combining the cases
    //in a smart way; kept them all separate for simplicity initially
    for (int k1 = Kmin(m, 0); k1<(K(0)+1); k1++){

      //Dominant species
      fg1 = calc_fg(k1, lam1(m), J, r1_sub, y1.row(m), miss.row(m));

      //-----------------------------------------------------------------------
      if((S == 2) & dep(1,0)){  //sp1 --> sp2

        lik_sp2 = 0;
        par2_m = lp2(m) + k1 * lp2_sp1(m);
        if(modOcc(1)){
          lik_sp2 = lik_subord_occ(Kmin(m, 1), par2_m, J, r2_sub, y2.row(m), miss.row(m));
        } else {
          lik_sp2 = lik_subord_abun(Kmin(m, 1), K(1), par2_m, J, r2_sub, y2.row(m), miss.row(m));
        }
        lik_m += fg1 * lik_sp2;

      //-----------------------------------------------------------------------
      } else if((S == 3) & (dep(1,0) & dep(2,1))){ // sp1 --> sp2 --> sp3

        lik_sp2 = 0;

        par2_m = exp(lp2(m) + k1 * lp2_sp1(m));

        for (int k2 = Kmin(m, 1); k2<(K(1)+1); k2++){
          fg2 = calc_fg(k2, par2_m, J, r2_sub, y2.row(m), miss.row(m));

          lik_sp3 = 0;
          par3_m = lp3(m) + k2 * lp3_sp2(m);
          if(modOcc(2)){
            lik_sp3 = lik_subord_occ(Kmin(m, 2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
          } else {
            lik_sp3 = lik_subord_abun(Kmin(m, 2), K(2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
          }
          lik_sp2 += fg2 * lik_sp3;
        }

        lik_m += fg1 * lik_sp2;

      //-----------------------------------------------------------------------
      } else if((S == 3) & (dep(2,0) & dep(2,1))){ // sp1 --> sp3 <-- sp2

        lik_sp2 = 0;

        for (int k2 = Kmin(m, 1); k2<(K(1)+1); k2++){

          fg2 = calc_fg(k2, lam2(m), J, r2_sub, y2.row(m), miss.row(m));

          par3_m = lp3(m) + k1 * lp3_sp1(m) + k2 * lp3_sp2(m);
          lik_sp3 = 0;
          if(modOcc(2)){
            lik_sp3 = lik_subord_occ(Kmin(m, 2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
          } else {
            lik_sp3 = lik_subord_abun(Kmin(m, 2), K(2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
          }
          lik_sp2 += fg2 * lik_sp3;
        }

        lik_m += fg1 * lik_sp2;

      //-----------------------------------------------------------------------
      } else if((S == 3) & (dep(1,0) & dep(2,0))){ // sp2 <-- sp1 --> sp3


        par2_m = lp2(m) + k1 * lp2_sp1(m);
        lik_sp2 = 0;
        if(modOcc(1)){
          lik_sp2 = lik_subord_occ(Kmin(m, 1), par2_m, J, r2_sub, y2.row(m), miss.row(m));
        } else {
          lik_sp2 = lik_subord_abun(Kmin(m, 1), K(1), par2_m, J, r2_sub, y2.row(m), miss.row(m));
        }

        par3_m = lp3(m) + k1 * lp3_sp1(m);
        lik_sp3 = 0;
        if(modOcc(2)){
          lik_sp3 = lik_subord_occ(Kmin(m, 2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
        } else {
          lik_sp3 = lik_subord_abun(Kmin(m, 2), K(2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
        }

        lik_m += fg1 * lik_sp2 * lik_sp3;
      } else if((S == 3) & (dep(1,0) & dep(2,0) & dep(2,1))){ // sp1 --> sp2 --> sp3 & sp1 --> sp3

        lik_sp2 = 0;

        par2_m = exp(lp2(m) + k1 * lp2_sp1(m));

        for (int k2 = Kmin(m, 1); k2<(K(1)+1); k2++){
          fg2 = calc_fg(k2, par2_m, J, r2_sub, y2.row(m), miss.row(m));

          lik_sp3 = 0;
          par3_m = lp3(m) + k1 * lp2_sp1(m) + k2 * lp3_sp2(m);
          if(modOcc(2)){
            lik_sp3 = lik_subord_occ(Kmin(m, 2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
          } else {
            lik_sp3 = lik_subord_abun(Kmin(m, 2), K(2), par3_m, J, r3_sub, y3.row(m), miss.row(m));
          }
          lik_sp2 += fg2 * lik_sp3;
        }

        lik_m += fg1 * lik_sp2;

      }
    }

    nll -= log(lik_m);
  }
  return nll;
}
