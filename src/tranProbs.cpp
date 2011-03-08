#include "tranProbs.h"

using namespace Rcpp ;

SEXP tranProbs( SEXP Nr, SEXP omegaR, SEXP gammaR, SEXP deltaR, SEXP dynamicsR)
{
    IntegerVector N(Nr);
    int n = N.size();
    double omega = as<double>(omegaR);
    double gamma = as<double>(gammaR);
    int delta = as<int>(deltaR);
    std::string dynamics = as<std::string>(dynamicsR) ;
    
    arma::mat bpsum = arma::zeros(n, n);
    for(int j=0; j<n; j++) {
        for(int k=0; k<n; k++) {
            int Nm = std::min(N[j], N[k]);
            IntegerVector cmin0 = seq(0, Nm);
            if(dynamics=="autoreg")
                bpsum(k, j) = sum(dbinom(cmin0, N[j], omega, false) *
                    dpois(N[k]-cmin0, gamma*N[j], false));
            else
                bpsum(k, j) = sum(dbinom(cmin0, N[j], omega, false) *
                    dpois(N[k]-cmin0, gamma, false));
            }
        }
    arma::mat cs1 = sum(bpsum, 0);
    arma::mat csm1 = arma::repmat(cs1, n, 1);
    bpsum = bpsum / csm1;
    if(delta > 1) {    
        for(int d=1; d<delta; d++) {
            bpsum *= bpsum;
            arma::mat cs = sum(bpsum, 0);
            arma::mat csm = arma::repmat(cs, n, 1);
            bpsum = bpsum / csm;
            }
        }
    return wrap(bpsum);
}

