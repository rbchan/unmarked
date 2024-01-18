#include <RcppArmadillo.h>
#include <float.h>
#include "detfuns.h"

// [[Rcpp::export]]
double nll_distsamp(Rcpp::IntegerMatrix y, Rcpp::NumericVector lam,
    Rcpp::NumericVector sig, double scale, Rcpp::NumericMatrix a,
    Rcpp::NumericMatrix u, Rcpp::NumericVector w, Rcpp::NumericVector db,
    std::string keyfun, std::string survey) {

  int R = y.nrow();   //y.n_rows;
  int J = y.ncol();   // y.n_cols;

  //  double mu = 0.0;
  double ll = 0.0;
  double lnmin = log(DBL_MIN);

  double f0 = 0.0;

  for(int i=0; i<R; i++) {
    if((survey=="line") && (keyfun=="halfnorm"))
      f0 = Rf_dnorm4(0.0, 0.0, sig[i], false);
    if((survey=="line") && (keyfun=="exp"))
      f0 = Rf_dexp(0.0, 1/sig[i], false);
    for(int j=0; j<J; j++) {
      double cp = 0.0;
      double lower = db[j];
      double upper = db[j+1];
      double result = 0.0;

      if(keyfun=="uniform") {
	cp = u(i,j);
      } else {
	if(survey=="point") {
	  if(keyfun=="halfnorm") {
	    result = sig[i]*sig[i]*(1-exp(-upper*upper/(2*sig[i]*sig[i])))-
	      sig[i]*sig[i]*(1-exp(-lower*lower/(2*sig[i]*sig[i])));
	  } else if(keyfun=="exp") {
      DetExp f(sig[i], 1);
      result = trap_rule(f, lower, upper);
	  } else if(keyfun=="hazard") {
      DetHaz f(sig[i], scale, 1);
      result = trap_rule(f, lower, upper);
	  }
	  cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	} else if(survey=="line") {
	  if(keyfun=="halfnorm") {
	    result = (Rf_pnorm5(upper, 0.0, sig[i], true, false) -
		      Rf_pnorm5(lower, 0.0, sig[i], true, false)) / f0;
	  } else if(keyfun=="exp") {
	    result = sig[i]*(1-exp(-upper/sig[i])) -
	      sig[i]*(1-exp(-lower/sig[i]));
	  } else if(keyfun=="hazard") {
      DetHaz f(sig[i], scale, 0);
      result = trap_rule(f, lower, upper);
	  }
	  cp = result / w[j] * u(i,j);
	}
      }
      ll += std::max(Rf_dpois(y(i,j), lam[i]*cp, true), lnmin);
    }
  }
  return -ll;
}
