#ifndef ANTL_VECTOR_FUNCTIONS_CPP
#define ANTL_VECTOR_FUNCTIONS_CPP



template<typename PType>
bool positive_func(PType dubby){
    if(dubby > 0) {return true;}
    else{ return false;}
};

template<typename PType>
bool is_close(const std::vector<PType> & vec1, const std::vector<PType> & vec2, const PType & maxdist){
    if(vec1.size() != vec2.size()){
      throw "vectors are not the same size";
    }
    PType absdif;

    for (int i = 0; i < vec1.size(); ++i ){
      sub(absdif, vec1[i], vec2[i]);
      abs(absdif, absdif);
      if(absdif > maxdist){
        return false;
      }
    }
    return true;
};

template<typename PType>
void compute_initial_s(const std::vector<PType> & alpha, const int kbound){
    int kprime;
    PType inf_norm = infinity_norm(alpha);
    std::vector<PType> s_term(alpha.size());
};


template<typename PType>
void log_to_valuation(std::vector<PType> &valuationvec, const std::vector<PType> & log_vec, const int r1){
  PType tempvar, new_coord = 0;
  valuationvec.clear();

  for(int i = 0; i < log_vec.size(); ++i ){
    if(i < r1){
      sub(new_coord, new_coord, log_vec[i]); //      new_coord -= log_vec[i];

      valuationvec.push_back(exp(log_vec[i]));
    }else{
      mul(tempvar, 2, log_vec[i]);
      sub(new_coord, new_coord, tempvar);      //new_coord -= 2*log_vec[i];
    }
  }
  if(log_vec.size()+1 > r1){
    div(new_coord, new_coord, 2);             //new_coord = new_coord/2;

  }
  valuationvec.push_back(exp(new_coord));

};






#endif // guard
