#include "detfuns.h"
#include <Rmath.h>

// Half-normal detection function for line-transects
void gxhn(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double sigma = v[0];
  for(int i=0; i<n; i++) {
    x[i] = exp(-x[i]*x[i] / (2*sigma*sigma));
  }
}

// Half-normal detection function for point-transects
void grhn(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double sigma = v[0];
  for(int i=0; i<n; i++) {
    x[i] = exp(-x[i]*x[i] / (2*sigma*sigma)) * x[i];
  }
}



// Negative exponential detection function for line-transects
void gxexp(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double rate = v[0];
  for(int i=0; i<n; i++) {
    x[i] = exp(-x[i] / rate);
  }
}

// Negative exponential detection function for point-transects
void grexp(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double rate = v[0];
  for(int i=0; i<n; i++) {
    x[i] = exp(-x[i] / rate) * x[i];
  }
}



// Hazard-rate detection function for line-transects
void gxhaz(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double shape = v[0];
  double scale = v[1];
  for(int i=0; i<n; i++) {
    x[i] = 1 - exp(-1*pow(x[i]/shape, -scale));
  }
}


// Hazard-rate detection function for point-transects
void grhaz(double *x, int n, void *ex) {
  double *v;
  v = (double*)ex; // cast to double pointer
  double shape = v[0];
  double scale = v[1];
  for(int i=0; i<n; i++) {
    x[i] = (1 - exp(-1*pow(x[i]/shape, -scale))) * x[i];
  }
}



