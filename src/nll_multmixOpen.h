#ifndef _unmarked_NLL_MULTMIXOPEN_H
#define _unmarked_NLL_MULTMIXOPEN_H

#include <RcppArmadillo.h>
#include <float.h>
#include "tranprobs.h"
#include "distr.h"
#include "pifun.h"
#include "utils.h"

RcppExport SEXP nll_multmixOpen( SEXP y_, SEXP yt_, SEXP Xlam_, SEXP Xgam_, SEXP Xom_,
    SEXP Xp_, SEXP Xiota_,
    SEXP beta_, SEXP beta_ind_,
    SEXP Xlam_offset_,
    SEXP Xgam_offset_, SEXP Xom_offset_, SEXP Xp_offset_, SEXP Xiota_offset_,
    SEXP ytna_, SEXP yna_, SEXP lk_, SEXP mixture_, 
    SEXP first_, SEXP last_, SEXP first1_,
    SEXP M_, SEXP T_, SEXP J_, SEXP R_, SEXP delta_, SEXP dynamics_,
    SEXP fix_, SEXP go_dims_, SEXP immigration_, 
    SEXP I_, SEXP I1_, SEXP Ib_, SEXP Ip_, SEXP pi_fun_,
    SEXP lfac_k_, SEXP kmyt_, SEXP lfac_kmyt_, SEXP fin_) ;

#endif
