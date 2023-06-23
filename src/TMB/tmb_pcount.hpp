#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

template<class Type>
Type lp_site_pcount(vector<Type> y, int mixture, Type lam, vector<Type> p, 
                    Type log_alpha, int K, int Kmin){
  
  Type alpha, var, f, g, out = 0.0;
  if(mixture == 2){
    alpha = exp(log_alpha);
    var = lam + pow(lam, 2) / alpha; //translate to TMB parameterization
  } else if(mixture == 3){
    alpha = invlogit(log_alpha);
  }
  
  for (int k=Kmin; k<(K+1); k++){
    if(mixture == 2){
      f = dnbinom2(Type(k), lam, var, false);
    } else if(mixture == 3){
      f = dzipois(Type(k), lam, alpha, false); 
    } else {
      f = dpois(Type(k), lam, false);
    }
    g = 0.0;
    for (int j=0; j<y.size(); j++){
      if(R_IsNA(asDouble(y(j)))) continue;
      g += dbinom(y(j), Type(k), p(j), true);
    }
    out += f * exp(g);
  }
  return log(out + DBL_MIN);
}

// name of function below **MUST** match filename
template <class Type>
Type tmb_pcount(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations
  DATA_INTEGER(K); //Max value of abundance for marginalization
  DATA_IVECTOR(Kmin); //Minimum obs at each site

  DATA_INTEGER(mixture); //Abundance distribution

  DATA_MATRIX(X_state); //lambda fixed effect design mat
  DATA_SPARSE_MATRIX(Z_state); //psi random effect design mat
  DATA_VECTOR(offset_state);
  DATA_INTEGER(n_group_vars_state); //# of grouping variables for lambda
  DATA_IVECTOR(n_grouplevels_state); //# of levels of each grouping variable

  DATA_MATRIX(X_det); //same thing but for p
  DATA_SPARSE_MATRIX(Z_det);
  DATA_VECTOR(offset_det);
  DATA_INTEGER(n_group_vars_det);
  DATA_IVECTOR(n_grouplevels_det);

  PARAMETER_VECTOR(beta_state); //Fixed effect params for lamda
  PARAMETER_VECTOR(b_state); //Random intercepts and/or slopes for lambda
  PARAMETER_VECTOR(lsigma_state); //Random effect variance(s) for lambda

  PARAMETER_VECTOR(beta_det); //Same thing but for det
  PARAMETER_VECTOR(b_det);
  PARAMETER_VECTOR(lsigma_det);

  PARAMETER_VECTOR(beta_scale);
  Type scale = 0;
  if(mixture > 1) scale = beta_scale(0);

  Type loglik = 0.0;
  int M = y.rows(); //# of sites
  int J = y.cols(); //# of observations per site

  //Construct lambda vector
  vector<Type> lam = X_state * beta_state + offset_state;
  lam = add_ranef(lam, loglik, b_state, Z_state, lsigma_state, 
                  n_group_vars_state, n_grouplevels_state);
  lam = exp(lam);

  //Construct p vector
  vector<Type> p = X_det * beta_det + offset_det;
  p = add_ranef(p, loglik, b_det, Z_det, lsigma_det, 
                n_group_vars_det, n_grouplevels_det);
  p = invlogit(p);
  
  //Likelihood
  for (int i=0; i<M; i++){
    int pstart = i * J;
    vector<Type> ysub = y.row(i);
    vector<Type> psub = p.segment(pstart, J);
    loglik -= lp_site_pcount(ysub, mixture, lam(i), psub,
                             scale, K, Kmin(i));
  }

  return loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
