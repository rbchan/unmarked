#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_occu(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations
  DATA_VECTOR(no_detect); //Indicator for if site had no detections

  DATA_INTEGER(link);

  DATA_MATRIX(X_state); //psi fixed effect design mat (state=psi)
  DATA_SPARSE_MATRIX(Z_state); //psi random effect design mat
  DATA_VECTOR(offset_state);
  DATA_INTEGER(n_group_vars_state); //# of grouping variables for psi
  DATA_IVECTOR(n_grouplevels_state); //# of levels of each grouping variable

  DATA_MATRIX(X_det); //same thing but for p
  DATA_SPARSE_MATRIX(Z_det);
  DATA_VECTOR(offset_det);
  DATA_INTEGER(n_group_vars_det);
  DATA_IVECTOR(n_grouplevels_det);

  PARAMETER_VECTOR(beta_state); //Fixed effect params for psi
  PARAMETER_VECTOR(b_state); //Random intercepts and/or slopes for psi
  PARAMETER_VECTOR(lsigma_state); //Random effect variance(s) for psi

  PARAMETER_VECTOR(beta_det); //Same thing but for det
  PARAMETER_VECTOR(b_det);
  PARAMETER_VECTOR(lsigma_det);

  Type loglik = 0.0;

  int M = y.rows(); //# of sites
  int J = y.cols(); //# of observations per site

  //Construct psi vector
  vector<Type> psi = X_state * beta_state + offset_state;
  psi = add_ranef(psi, loglik, b_state, Z_state, lsigma_state, 
                  n_group_vars_state, n_grouplevels_state);
  if(link == 1){
    psi = cloglog(psi);
  } else {
    psi = invlogit(psi);
  }

  //Construct p vector
  vector<Type> p = X_det * beta_det + offset_det;
  p = add_ranef(p, loglik, b_det, Z_det, lsigma_det, 
                n_group_vars_det, n_grouplevels_det);
  p = invlogit(p);

  //Standard occupancy likelihood calculation
  for (int i=0; i<M; i++){
    int pind = i * J;
    Type cp = 1.0;
    for (int j=0; j<J; j++){
      if(R_IsNA(asDouble(y(i,j)))){
        pind += 1;
        continue;
      }
      cp *= pow(p(pind), y(i,j)) * pow(1-p(pind), 1-y(i,j));
      pind += 1;
    }
    loglik -= log(psi(i) * cp + (1-psi(i)) * no_detect(i));
  }

  return loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

