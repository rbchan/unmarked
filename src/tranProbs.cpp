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
            if(dynamics=="autoreg") {
                for(int c=0; c<Nm; c++) 
                    bpsum(j, k) += sum(Rf_dbinom(c, N[j], omega, false) *
                        Rf_dpois(N[k]-c, gamma*N[j], false));
                }
            else {
                for(int c=0; c<Nm; c++)
                    bpsum(j, k) += sum(Rf_dbinom(c, N[j], omega, false) *
                        Rf_dpois(N[k]-c, gamma, false));
                }
            }
        }
    if(delta > 1) {    
        for(int d=1; d<delta; d++) {
            bpsum *= bpsum;
//            arma::mat cs = sum(bpsum, 0);
//            arma::mat csm = arma::repmat(cs, n, 1);
//            bpsum = bpsum / csm;
            }
        }
    return wrap(bpsum);
}

