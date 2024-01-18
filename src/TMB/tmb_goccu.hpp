#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// Adapted from Stan code by Maxwell B. Joseph,
// https://discourse.mc-stan.org/t/divergent-transition-every-iteration-multi-scale-occupancy-model/13739/5

// name of function below **MUST** match filename
template<class Type>
Type tmb_goccu(objective_function<Type>* obj) {
  //Describe input data
  DATA_MATRIX(y); //observations
  DATA_INTEGER(T);

  DATA_INTEGER(link);

  DATA_MATRIX(Xpsi);
  DATA_MATRIX(Xphi);
  DATA_MATRIX(Xp);

  DATA_INTEGER(n_possible);
  DATA_MATRIX(alpha_potential);
  DATA_VECTOR(known_present);
  DATA_MATRIX(known_available);
  DATA_MATRIX(missing_session);

  PARAMETER_VECTOR(beta_psi);
  PARAMETER_VECTOR(beta_phi);
  PARAMETER_VECTOR(beta_p);

  vector<Type> psi = Xpsi * beta_psi;
  vector<Type> phi = Xphi * beta_phi;
  vector<Type> p = Xp * beta_p;

  if(link == 1){
    //psi = cloglog(psi);
  } else {
    psi = invlogit(psi);
  }

  phi = invlogit(phi);
  p = invlogit(p);

  Type loglik = 0.0;

  int M = y.rows();
  int J = y.cols() / T;

  Type obs_lp;
  Type poss_lp;
  Type exp_poss_lp;

  int ystart;
  vector<Type> ysite;
  vector<Type> psite;
  vector<Type> ysub;
  vector<Type> psub;

  int tstart;
  int pstart;

  for (int i=0; i<M; i++){

    tstart = i*T;
    pstart = i*(T*J);

    ysite = y.row(i);
    psite = p.segment(pstart, T*J);

    if(known_present(i) == 1){
      for (int t=0; t<T; t++){

        if(missing_session(i, t) == 1) continue;

        ystart = t*J;

        psub = psite.segment(ystart, J);
        ysub = ysite.segment(ystart, J);

        obs_lp = log(phi(tstart+t));
        for (int j=0; j<J; j++){
            if(R_IsNA(asDouble(ysub(j)))) continue;
            obs_lp += dbinom(ysub(j), Type(1), psub(j), true);
        }
        if(known_available(i, t) == 1){
          loglik += obs_lp;
        } else {
          loglik += log(exp(obs_lp) + exp(log(1-phi(tstart+t))));
        }
      }
      loglik += log(psi(i));
    } else {

      exp_poss_lp = 0.0;

      for (int k=0; k<n_possible; k++){
        poss_lp = log(psi(i));

        for (int t=0; t<T; t++){

          if(missing_session(i, t) == 1) continue;

          ystart = t*J;

          psub = psite.segment(ystart, J);
          ysub = ysite.segment(ystart, J);

          if(alpha_potential(k, t) == 0){
            poss_lp += log(1 - phi(tstart+t));
          } else {
            poss_lp += log(phi(tstart+t));
            for (int j=0; j<J; j++){
              if(R_IsNA(asDouble(ysub(j)))) continue;
              poss_lp += dbinom(ysub(j), Type(1), psub(j), true);
            }
          }
        }
        exp_poss_lp += exp(poss_lp);
      }
      exp_poss_lp += exp(log(1-psi(i)));

    loglik += log(exp_poss_lp);
    }
  }

  return -loglik;

}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
