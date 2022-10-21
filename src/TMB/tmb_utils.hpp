template<class Type>
vector<Type> cloglog(vector<Type> inp) {
  int sz = inp.size();
  vector<Type> out(sz);
  for (int i=0; i<sz; i++){
    out(i) = 1 - exp(-exp(inp(i)));
  }
  return out;
}

template<class Type>
vector<Type> add_ranef(vector<Type> par, Type& loglik, 
                vector<Type> b, Eigen::SparseMatrix<Type> Z,
                vector<Type> lsigma, int n_group_vars, vector<int> n_grouplevels) {
  
  if(n_group_vars == 0) return par;
  vector<Type> sigma = exp(lsigma);
  int idx = 0;
  for (int i=0; i<n_group_vars; i++){
    for (int j=0; j<n_grouplevels(i); j++){
      loglik -= dnorm(b(idx), Type(0.0), sigma(i), true);
      idx += 1;
    }
  }
  par += Z * b;
  return par;
}
