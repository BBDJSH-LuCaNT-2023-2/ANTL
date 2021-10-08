#ifndef ANTL_VECTOR_FUNCTIONS_CPP
#define ANTL_VECTOR_FUNCTIONS_CPP

template<typename PType>
bool is_close(std::vector<PType> & vec1, std::vector<PType> & vec2, PType & maxdist){
    if(vec1.size() != vec2.size()){
      throw "vectors are not the same size";
    }
    for (int i = 0; i < vec1.size(); ++i ){
      if(abs(vec1[i] -vec2[i] ) > maxdist){
        return false;
      }
    }
    return true;
};





#endif // guard
