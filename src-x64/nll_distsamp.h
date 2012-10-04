#ifndef _UNMARKED_NLL_DISTSAMP_H
#define _UNMARKED_NLL_DISTSAMP_H

#include "R_ext/Applic.h"

#include <Rcpp.h>


RcppExport SEXP nll_distsamp( SEXP y_, SEXP lam_, SEXP sig_, SEXP scale_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_, SEXP survey_, SEXP reltol_ ) ;

#endif
