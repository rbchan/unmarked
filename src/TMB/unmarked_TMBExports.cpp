#define TMB_LIB_INIT R_init_unmarked_TMBExports
#include <float.h>
#include <TMB.hpp>
#include "tmb_utils.hpp"
#include "tmb_pifun.hpp"
#include "tmb_keyfun.hpp"
#include "tmb_occu.hpp"
#include "tmb_pcount.hpp"
#include "tmb_multinomPois.hpp"
#include "tmb_distsamp.hpp"
#include "tmb_gdistremoval.hpp"
#include "tmb_IDS.hpp"
#include "tmb_goccu.hpp"

template<class Type>
Type objective_function<Type>::operator() () {
  DATA_STRING(model);
  if(model == "tmb_occu") {
    return tmb_occu(this);
  } else if(model == "tmb_pcount") {
    return tmb_pcount(this);
  } else if(model == "tmb_multinomPois"){
    return tmb_multinomPois(this);
  } else if(model == "tmb_distsamp"){
    return tmb_distsamp(this);
  } else if(model == "tmb_gdistremoval"){
    return tmb_gdistremoval(this);
  } else if(model == "tmb_IDS"){
    return tmb_IDS(this);
  } else if(model == "tmb_goccu"){
    return tmb_goccu(this);
  } else {
    error("Unknown model.");
  }
  return 0;
}
