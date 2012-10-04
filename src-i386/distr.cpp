
#include <Rmath.h>
#include "distr.h"

// Zero-inflated Poisson
double dzip(int x, double lambda, double psi) {
  double den=0.0;
  if(x==0)
    den = psi + (1-psi)*exp(-lambda);
  else if(x>0)
    den = (1-psi)*Rf_dpois(x, lambda, false);
  return den;
}

