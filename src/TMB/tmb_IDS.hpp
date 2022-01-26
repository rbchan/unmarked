#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_IDS(objective_function<Type>* obj) {
 DATA_MATRIX(pind);
  DATA_VECTOR(lam_adjust);

  DATA_MATRIX(y_hds);
  DATA_MATRIX(X_hds);
  DATA_MATRIX(V_hds);
  DATA_INTEGER(key_hds);
  DATA_VECTOR(db_hds);
  DATA_VECTOR(a_hds);
  DATA_VECTOR(w_hds);
  DATA_VECTOR(u_hds);

  DATA_MATRIX(y_pc);
  DATA_MATRIX(X_pc);
  DATA_MATRIX(V_pc);
  DATA_INTEGER(key_pc);
  DATA_VECTOR(db_pc);
  DATA_VECTOR(a_pc);
  DATA_VECTOR(w_pc);
  DATA_VECTOR(u_pc);

  DATA_MATRIX(y_oc);
  DATA_MATRIX(X_oc);
  DATA_MATRIX(V_oc);
  DATA_INTEGER(key_oc);
  DATA_VECTOR(db_oc);
  DATA_VECTOR(a_oc);
  DATA_VECTOR(w_oc);
  DATA_VECTOR(u_oc);
  DATA_INTEGER(K_oc);
  DATA_IVECTOR(Kmin_oc);

  DATA_VECTOR(durationDS);
  DATA_VECTOR(durationPC);
  DATA_VECTOR(durationOC);
  DATA_MATRIX(Xa_hds);
  DATA_MATRIX(Xa_pc);
  DATA_MATRIX(Xa_oc);

  PARAMETER_VECTOR(beta_lam);
  PARAMETER_VECTOR(beta_hds);
  PARAMETER_VECTOR(beta_pc);
  PARAMETER_VECTOR(beta_oc);
  PARAMETER_VECTOR(beta_avail);

  int survey = 1; // Only point counts supported

  Type loglik = 0;

  // HDS
  int M = y_hds.rows();
  int J = y_hds.cols();

  vector<Type> lam_hds = X_hds * beta_lam;
  lam_hds = exp(lam_hds);
  lam_hds = lam_hds * lam_adjust(0);

  vector<Type> sigma_hds = V_hds * beta_hds;
  sigma_hds = exp(sigma_hds);

  vector<Type> p_avail(M);
  p_avail.setOnes();
  if(beta_avail.size() > 0){
    vector<Type> phi = Xa_hds * beta_avail;
    phi = exp(phi);
    p_avail = 1 - exp(-1 * durationDS.array() * phi.array());
  }

  for (int i=0; i<M; i++){

    vector<Type> cp = distance_prob(key_hds, sigma_hds(i), Type(0), survey,
                                    db_hds, w_hds, a_hds, u_hds);

    for (int j=0; j<J; j++){
      loglik -= dpois(y_hds(i,j), lam_hds(i) * cp(j) * p_avail(i), true);
    }
  }

  //printf("hds done");

  // Point count
  M = y_pc.rows();

  if(M > 0){

  vector<Type> lam_pc = X_pc * beta_lam;
  lam_pc = exp(lam_pc);
  lam_pc = lam_pc * lam_adjust(1);

  //vector<Type> sigma_pc = V_pc * beta_pc;
  vector<Type> sigma_pc;
  if(beta_pc.size() > 0){
    sigma_pc = V_pc * beta_pc;
  } else{
    sigma_pc = V_pc * beta_hds;
  }
  sigma_pc = exp(sigma_pc);

  vector<Type> p_avail_pc(M);
  p_avail_pc.setOnes();
  if(beta_avail.size() > 0){
    vector<Type> phi = Xa_pc * beta_avail;
    phi = exp(phi);
    p_avail_pc = 1 - exp(-1 * durationPC.array() * phi.array());
  }

  for (int i=0; i<M; i++){
    vector<Type> cp = distance_prob(key_pc, sigma_pc(i), Type(0), survey,
                                    db_pc, w_pc, a_pc, u_pc);
    loglik -= dpois(y_pc(i,0), lam_pc(i) * cp(0) * p_avail_pc(i), true);
  }

  }

  //printf("pc done");

  // R-N occupancy
  M = y_oc.rows();

  if(M > 0){

  vector<Type> lam_oc = X_oc * beta_lam;
  lam_oc = exp(lam_oc);
  lam_oc = lam_oc * lam_adjust(2);

  //vector<Type> sigma_oc = V_oc * beta_oc;
  vector<Type> sigma_oc;
  if(beta_oc.size() > 0){
    sigma_oc = V_oc * beta_oc;
  } else{
    sigma_oc = V_oc * beta_hds;
  }
  sigma_oc = exp(sigma_oc);

  vector<Type> p_avail_oc(M);
  p_avail_oc.setOnes();
  if(beta_avail.size() > 0){
    vector<Type> phi = Xa_oc * beta_avail;
    phi = exp(phi);
    p_avail_oc = 1 - exp(-1 * durationOC.array() * phi.array());
  }

  Type f;
  Type g;
  Type p;
  Type site_lp;

  for (int i=0; i<M; i++){

    vector<Type> q = 1 - distance_prob(key_oc, sigma_oc(i), Type(0), survey,
                                    db_oc, w_oc, a_oc, u_oc) * p_avail_oc(i);
    //vector<Type> q = 1 - distance_prob(key_oc, sigma_oc(i), Type(0), survey,
    //                                db_oc, w_oc, a_oc, u_oc);

    site_lp = 0.0;
    for (int k=Kmin_oc(i); k<(K_oc+1); k++){
      f = dpois(Type(k), lam_oc(i), false);
      p = 1 - pow(q(0), Type(k));
      g = dbinom(y_oc(i,0), Type(1), p, false);
      site_lp += f * g;
    }
    loglik -= log(site_lp);

  }

  }

  return loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
