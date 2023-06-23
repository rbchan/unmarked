//Base class for functions to integrate
template<class Type>
class IntFunc {
  public:
    IntFunc() {}
    
    virtual Type operator()(const Type& x) const {
      return(x);
    }
};

//Negative exponential function
template<class Type>
class DetExp: public IntFunc<Type> {
  private:
    Type rate;
    int point;
  public:
    DetExp(Type rate_, int point_) :
      rate(rate_), point(point_) {}

    Type operator()(const Type& x) const {
      Type pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return( exp( -x/rate) * pd_adjust);
    }
};

//Hazard function
template<class Type>
class DetHaz: public IntFunc<Type> {
  private:
    Type shape;
    Type scale;
    int point;
  public:
    DetHaz(Type shape_, Type scale_, int point_) : 
      shape(shape_), scale(scale_), point(point_) {}

    Type operator()(const Type& x) const {
      Type pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return( (1-exp(-1 * pow(x/shape, -scale))) * pd_adjust);
    }
};

//Integrate with trapezoidal rule
template<class Type>
Type trap_rule(IntFunc<Type> &f, Type a, Type b){
  
  int n = 100;
  Type h = (b - a) / n;
  
  Type int_sum = 0;
  for(int i=1; i<n; i++){
    int_sum += f(a + i*h);
  }

  return( h/2 * (f(a) + 2*int_sum + f(b)) );
}

template<class Type>
vector<Type> key_halfnorm(Type sigma, int survey_type, vector<Type> db, 
                          vector<Type> w, vector<Type> a){

  int J = db.size() - 1;
  vector<Type> p(J);
  Type p1, p2;

  if(survey_type == 0){
    Type f0 = 2 * dnorm(Type(0.0), Type(0.0), sigma, false);
    for (int j=0; j<J; j++){
      p1 = pnorm(db(j+1), Type(0.0), sigma);
      p2 = pnorm(db(j), Type(0.0), sigma);
      p(j) = 2 * (p1 - p2) / f0 / w(j);
    }

  } else if(survey_type == 1){
    Type s2 = pow(sigma, 2);
    for (int j=0; j<J; j++){
      p1 = 1 - exp(-pow(db(j+1),2) / (2 * s2));
      p2 = 1 - exp(-pow(db(j),2) / (2 * s2)); 
      p(j) = (s2 * p1 - s2 * p2) * 2 * M_PI / a(j);
    }
  }
  return(p);
}

template<class Type>
vector<Type> key_exp(Type rate, int survey_type, vector<Type> db, 
                     vector<Type> w, vector<Type> a){

  int J = db.size() - 1;
  vector<Type> p(J);

  if(survey_type == 0){
    for(int j=0; j<J; j++){
      Type int_ = rate*(1 - exp(-db(j+1)/rate)) - rate*(1-exp(-db(j)/rate));
      p(j) = int_ / w(j);
    }

  } else if(survey_type == 1){
    DetExp<Type> f(rate, 1);
    for(int j=0; j<J; j++){
      Type int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

template<class Type>
vector<Type> key_hazard(Type shape, Type scale, int survey_type, 
                        vector<Type> db, vector<Type> w, vector<Type> a){

  int J = db.size() - 1;
  vector<Type> p(J);
  DetHaz<Type> f(shape, scale, survey_type);

  if(survey_type == 0){
    for(int j=0; j<J; j++){
      Type int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ / w(j);
    }

  } else if(survey_type == 1){
    for(int j=0; j<J; j++){
      Type int_ = trap_rule(f, db(j), db(j+1));
      p(j) = int_ * 2 * M_PI / a(j);
    }
  }
  return(p);
}

template<class Type>
vector<Type> distance_prob(int keyfun_type, Type param1, Type param2, 
                           int survey_type, vector<Type> db,
                           vector<Type> w, vector<Type> a, vector<Type> u){
  
  int J = db.size() - 1;
  vector<Type> p(J);
  if(keyfun_type == 0){ // uniform
    p.setOnes(); 
  } else if (keyfun_type == 1){ // half-normal
    //param1 is sigma
    p = key_halfnorm(param1, survey_type, db, w, a);
  } else if (keyfun_type == 2){ // exponential
    //param1 is rate
    p = key_exp(param1, survey_type, db, w, a);
  } else if (keyfun_type == 3){ // hazard
    //param1 is shape, param2 is scale
    p = key_hazard(param1, param2, survey_type, db, w, a);
  } else{
    throw(std::invalid_argument("invalid keyfun"));
  }

  p = p.array() * u.array();

  return(p);

}
