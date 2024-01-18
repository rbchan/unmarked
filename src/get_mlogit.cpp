#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat get_mlogit(arma::mat lp_mat, std::string type, int S, arma::mat guide){

  int R = lp_mat.n_rows;
  int C = lp_mat.n_cols;

  if(type == "psi"){

    mat out = exp(lp_mat);
    for(int r=0; r<R; r++){
      double row_sum = sum(out.row(r)) + 1;
      for(int c=0; c<C; c++){
        out(r,c) = out(r,c) / row_sum;
      }
    }
    return out;

  } else if(type == "phi"){

    mat out(R,C);

    for(int r=0; r<R; r++){
      int ix = 0;
      rowvec p_row = exp(lp_mat.row(r));
      for(int s=0; s<S; s++){
        rowvec sub = p_row.subvec(ix,(ix+S-2));
        out(r,span(ix,(ix+S-2))) = sub / (sum(sub)+1);
        ix += (S-1);
      }
    }
    return out;
  } else if(type == "det"){

    mat out(R,C);
    for(int r=0; r<R; r++){
      mat sdp = zeros(S,S);
      sdp.col(0) = ones(S);
      for(int c=0; c<C; c++){
        sdp( guide(c,0), guide(c,1)) = exp(lp_mat(r,c));
      }
      for(int s=0; s<S; s++){
        double row_sum = sum(sdp.row(s));
        for (int j=0; j<S; j++){
          sdp(s,j) = sdp(s,j) / row_sum;
        }
      }
      for(int c=0; c<C; c++){
        out(r,c) = sdp( guide(c,0), guide(c,1) );
      }
    }
    return out;
  } else {
    stop("type must be 'psi','phi',or 'det'");
  }
}
