#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_distsamp(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations

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
  
  DATA_INTEGER(survey_type);
  DATA_INTEGER(keyfun_type);

  DATA_VECTOR(A); // Area
  DATA_VECTOR(db); // distance breaks
  DATA_MATRIX(a);
  DATA_VECTOR(w);
  DATA_MATRIX(u);

  PARAMETER_VECTOR(beta_state); //Fixed effect params for lambda
  PARAMETER_VECTOR(b_state); //Random intercepts and/or slopes for lambda
  PARAMETER_VECTOR(lsigma_state); //Random effect variance(s) for lambda

  PARAMETER_VECTOR(beta_det); //Same thing but for det
  PARAMETER_VECTOR(b_det);
  PARAMETER_VECTOR(lsigma_det);
 
  PARAMETER_VECTOR(beta_scale); //Trick here: this is 0-length array if keyfun != hazard
  Type scale = 0; // If not hazard  this is ignored later 
  if(keyfun_type == 3) scale = exp(beta_scale(0)); // If hazard

  Type loglik = 0.0;

  int M = y.rows(); // # of sites
  int J = y.cols(); // # of distance categories per site

  //Construct lambda vector
  vector<Type> lam = X_state * beta_state + offset_state;
  lam = add_ranef(lam, loglik, b_state, Z_state, lsigma_state, 
                  n_group_vars_state, n_grouplevels_state);
  lam = exp(lam);
  lam = lam.array() * A.array();

  //Construct distance parameter (sigma, rate, etc.) vector
  vector<Type> dp(M);
  if(keyfun_type > 0){ // If keyfun is not uniform
    dp = X_det * beta_det + offset_det;
    dp = add_ranef(dp, loglik, b_det, Z_det, lsigma_det, 
                   n_group_vars_det, n_grouplevels_det);
    dp = exp(dp);
  }

  //Likelihood
  for (int i=0; i<M; i++){
    //Not sure if defining this inside loop is necessary for parallel
    vector<Type> asub = a.row(i);
    vector<Type> usub = u.row(i);
    vector<Type> cp = distance_prob(keyfun_type, dp(i), scale, survey_type, db, 
                                    w, asub, usub); 
    vector<Type> ysub = y.row(i);
    Type site_lp = 0;
    
    for (int j=0; j<J; j++){
      if(R_IsNA(asDouble(ysub(j)))) goto endsite; //If any NAs found skip site
      site_lp += dpois(ysub(j), lam(i) * cp(j), true);
    }
    loglik -= site_lp;
    endsite: ;
  }

  return loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
