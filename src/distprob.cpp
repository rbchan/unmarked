#include "distprob.h"

using namespace Rcpp;
using namespace arma;

vec p_halfnorm(const double& sigma, const std::string& survey,
               const vec& db, const vec& w, const rowvec& a){

  int J = db.size() - 1;
  vec p(J);
  double p1, p2;

  if(survey == "line"){
    double f0 = 2 * R::dnorm(0.0, 0.0, sigma, 0);
    for (int j=0; j<J; j++){
      p1 = R::pnorm(db(j+1), 0.0, sigma, 1, 0);
      p2 = R::pnorm(db(j), 0.0, sigma, 1, 0);
      p(j) = 2 * (p1 - p2) / f0 / w(j);
    }

  } else if(survey == "point"){
    double s2 = pow(sigma, 2);
    for (int j=0; j<J; j++){
      p1 = 1 - exp(-pow(db(j+1),2) / (2 * s2));
      p2 = 1 - exp(-pow(db(j),2) / (2 * s2));
      p(j) = (s2 * p1 - s2 * p2) * 2 * M_PI / a(j);
    }
  }
  return(p);
}

vec p_exp(const double& rate, const std::string& survey, const vec& db,
          const vec& w, const rowvec& a, double& rel_tol){

  int J = db.size() - 1;
  vec p(J);

  if(survey == "line"){
    for(int j=0; j<J; j++){
      double int_ = rate*(1 - exp(-db(j+1)/rate)) - rate*(1-exp(-db(j)/rate));
      p(j) = int_ / w(j);
    }

  } else if(survey == "point"){
    DetExp f(rate, 1);
    for(int j=0; j<J; j++){
      double int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

vec p_hazard(const double& shape, const double& scale, const std::string& survey,
          const vec& db, const vec& w, const rowvec& a, double& rel_tol){

  int J = db.size() - 1;
  vec p(J);

  if(survey == "line"){
    DetHaz f(shape, scale, 0);
    for(int j=0; j<J; j++){
      double int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ / w(j);
    }

  } else if(survey == "point"){
    DetHaz f(shape, scale, 1);
    for(int j=0; j<J; j++){
      double int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

vec distprob(const std::string& keyfun, const double param1,
             const double param2, const std::string& survey, const vec& db,
             const vec& w, const rowvec& a){

  int J = db.size() - 1;
  double rel_tol = 0.0; //might use in future
  vec p(J);
  if(keyfun == "uniform"){
    p = ones(J);
  } else if (keyfun == "halfnorm"){
    //param1 is sigma
    p = p_halfnorm(param1, survey, db, w, a);
  } else if (keyfun == "exp"){
    //param1 is rate
    p = p_exp(param1, survey, db, w, a, rel_tol);
  } else if (keyfun == "hazard"){
    //param1 is shape, param2 is scale
    p = p_hazard(param1, param2, survey, db, w, a, rel_tol);
  } else{
    stop("invalid keyfun");
  }

  return(p);

}
