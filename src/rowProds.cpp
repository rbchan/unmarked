#include "rowProds.h"

using namespace Rcpp ;

SEXP rowProds( SEXP m )
{
    arma::mat M = as<arma::mat>(m);
    arma::colvec p = prod(M, 1);
    typedef std::vector<double> stdvec;
    stdvec out = arma::conv_to<stdvec>::from(p);
    return wrap(out);
}

