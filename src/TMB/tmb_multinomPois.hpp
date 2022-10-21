#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_multinomPois(objective_function<Type>* obj) {
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
  
  DATA_INTEGER(pifun_type);

  PARAMETER_VECTOR(beta_state); //Fixed effect params for lamda
  PARAMETER_VECTOR(b_state); //Random intercepts and/or slopes for lambda
  PARAMETER_VECTOR(lsigma_state); //Random effect variance(s) for lambda

  PARAMETER_VECTOR(beta_det); //Same thing but for det
  PARAMETER_VECTOR(b_det);
  PARAMETER_VECTOR(lsigma_det);

  //Define the log likelihood so that it can be calculated in parallel over sites
  Type loglik = 0.0;

  int M = y.rows(); // # of sites
  int J = y.cols(); // # of obs per site
  int R = X_det.rows() / M; // # of detection probs per site

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
    //Not sure if defining these inside loop is necessary for parallel
    int pstart = i * R;
    vector<Type> ysub = y.row(i);
    vector<Type> psub = p.segment(pstart, R);
    vector<Type> pi_lam = pifun(psub, pifun_type) * lam(i);
    
    for (int j=0; j<J; j++){
      if(R_IsNA(asDouble(ysub(j)))) continue;
      loglik -= dpois(ysub(j), pi_lam(j), true);
    }
  }

  return loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
