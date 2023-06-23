#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_gdistremoval(objective_function<Type>* obj) {
  //Describe input data
  DATA_VECTOR(y_dist); //observations
  DATA_VECTOR(y_rem);
  DATA_MATRIX(y_sum);
  DATA_INTEGER(mixture);
  DATA_INTEGER(K);
  DATA_IVECTOR(Kmin);
  DATA_INTEGER(T);

  DATA_MATRIX(X_lambda); //lambda fixed effect design mat
  DATA_SPARSE_MATRIX(Z_lambda); //psi random effect design mat
  DATA_INTEGER(n_group_vars_lambda); //# of grouping variables for lambda
  DATA_IVECTOR(n_grouplevels_lambda); //# of levels of each grouping variable
  
  DATA_MATRIX(X_phi);
  DATA_SPARSE_MATRIX(Z_phi);
  DATA_INTEGER(n_group_vars_phi);
  DATA_IVECTOR(n_grouplevels_phi);

  DATA_MATRIX(X_dist);
  DATA_SPARSE_MATRIX(Z_dist);
  DATA_INTEGER(n_group_vars_dist);
  DATA_IVECTOR(n_grouplevels_dist);

  DATA_MATRIX(X_rem);
  DATA_SPARSE_MATRIX(Z_rem);
  DATA_INTEGER(n_group_vars_rem);
  DATA_IVECTOR(n_grouplevels_rem);
  
  DATA_INTEGER(keyfun_type);

  DATA_VECTOR(A); // Area
  DATA_VECTOR(db); // distance breaks
  DATA_MATRIX(a);
  DATA_VECTOR(w);
  DATA_MATRIX(u);

  DATA_IVECTOR(per_len); // Length of removal periods

  PARAMETER_VECTOR(beta_lambda); //Fixed effect params for lambda
  PARAMETER_VECTOR(b_lambda); //Random intercepts and/or slopes for lambda
  PARAMETER_VECTOR(lsigma_lambda); //Random effect variance(s) for lambda
  
  PARAMETER_VECTOR(beta_alpha); //Only used if NB or ZIP
  Type log_alpha = 0;
  if(mixture > 1) log_alpha = beta_alpha(0);

  PARAMETER_VECTOR(beta_phi);
  PARAMETER_VECTOR(b_phi);
  PARAMETER_VECTOR(lsigma_phi);

  PARAMETER_VECTOR(beta_dist); //Same thing but for det
  PARAMETER_VECTOR(b_dist);
  PARAMETER_VECTOR(lsigma_dist);
 
  PARAMETER_VECTOR(beta_scale); //Trick here: this is 0-length array if keyfun != hazard
  Type scale = 0; // If not hazard  this is ignored later 
  if(keyfun_type == 3) scale = exp(beta_scale(0)); // If hazard

  PARAMETER_VECTOR(beta_rem); //Same thing but for det
  PARAMETER_VECTOR(b_rem);
  PARAMETER_VECTOR(lsigma_rem);
  
  Type loglik = 0.0;
  int M = X_lambda.rows(); // # of sites
  int Rdist = y_dist.size() / M;
  int Jdist = Rdist / T;
  int Rrem = y_rem.size() / M;
  int Jrem = Rrem / T;

  //Construct lambda vector
  vector<Type> lam = X_lambda * beta_lambda;
  lam = add_ranef(lam, loglik, b_lambda, Z_lambda, lsigma_lambda, 
                  n_group_vars_lambda, n_grouplevels_lambda);
  lam = exp(lam);
  lam = lam.array() * A.array();
  
  //Construct availability (phi) vector
  vector<Type> phi(M*T);
  phi.setOnes();
  if(T > 1){
    phi = X_phi * beta_phi;
    phi = add_ranef(phi, loglik, b_phi, Z_phi, lsigma_phi,
                    n_group_vars_phi, n_grouplevels_phi);
    phi = invlogit(phi);
  }

  //Construct distance parameter (sigma, rate, etc.) vector
  vector<Type> dp(M*T);
  if(keyfun_type > 0){ // If keyfun is not uniform
    dp = X_dist * beta_dist;
    dp = add_ranef(dp, loglik, b_dist, Z_dist, lsigma_dist, 
                   n_group_vars_dist, n_grouplevels_dist);
    dp = exp(dp);
  }


  //Construct removal parameter vector
  vector<Type> rp(M*Rrem);
  rp = X_rem * beta_rem;
  rp = add_ranef(rp, loglik, b_rem, Z_rem, lsigma_rem,
                 n_group_vars_rem, n_grouplevels_rem);
  rp = invlogit(rp);

  //Likelihood
  for (int i=0; i<M; i++){
    Type site_lp = 0;
    
    // Calculate 2nd abundance parameter
    Type alpha = 0;
    Type var;
    if(mixture == 2){
      alpha = exp(log_alpha);
      var = lam(i) + pow(lam(i), 2) / alpha; // TMB parameterization
    } else if(mixture == 3){
      alpha = invlogit(alpha);
    }

    vector<Type> f(K+1);
    f.setZero();
    // Iterate over possible true abundance values
    for (int k=Kmin(i); k<(K+1); k++){
      if(mixture == 2){
        f(k) = dnbinom2(Type(k), lam(i), var, false);
      } else if(mixture == 3){
        f(k) = dzipois(Type(k), lam(i), alpha, false);
      } else {
        f(k) = dpois(Type(k), lam(i), false);
      }
    }

    int t_ind = i * T;
    int yd_ind = i * T * Jdist;
    int yr_ind = i * T * Jrem;

    vector<Type> yd_sub(Jdist);
    vector<Type> yr_sub(Jrem);
    vector<Type> rp_sub(Jrem);
    vector<Type> cpd(Jdist);
    vector<Type> cpr(Jrem);
    Type pdist;
    Type prem;
    Type pall;
    vector<Type> g(K+1);
    g.setOnes();
    vector<Type> fg(K+1);

    vector<Type> asub = a.row(i);
    vector<Type> usub = u.row(i);

    for (int t=0; t<T; t++){
      yd_sub = y_dist.segment(yd_ind, Jdist);
      yr_sub = y_rem.segment(yr_ind, Jrem);
      rp_sub = rp.segment(yr_ind, Jrem);

      //If any NAs found skip site
      for (int j=0; j<Jdist; j++){
        if(R_IsNA(asDouble(yd_sub(j)))) goto endper; 
      }
      for (int j=0; j<Jrem; j++){
        if(R_IsNA(asDouble(yr_sub(j)))) goto endper;
      }

      cpd = distance_prob(keyfun_type, dp(t_ind), scale, 1, db, w, asub, usub);
      pdist = sum(cpd);
      cpd = cpd/pdist;
      site_lp += dmultinom(yd_sub, cpd, true);

      cpr = pifun_removal(rp_sub, per_len);
      prem = sum(cpr);
      cpr = cpr/sum(cpr);
      site_lp += dmultinom(yr_sub, cpr, true);

      pall = pdist * prem * phi(t_ind);
      
      for (int k=Kmin(i); k<(K+1); k++){
        g(k) *= dbinom(y_sum(i,t), Type(k), pall, false);
      }
      
      endper: ;
      t_ind += 1;
      yd_ind += Jdist;
      yr_ind += Jrem;
    }
    
    fg = f.array() * g.array();
    site_lp += log(sum(fg));
    
    loglik -= site_lp;
  }

  return loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
