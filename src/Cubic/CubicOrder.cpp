#ifndef ANTL_CUBIC_ORDER_CPP
#define ANTL_CUBIC_ORDER_CPP

#include "../../include/ANTL/Cubic/CubicOrder.hpp"


// constructor definitions
template <typename Type, typename PType>
CubicOrder<Type, PType> :: CubicOrder(polynomial<Type> const &poly) {
    defining_IBCF = poly;



}




// definitions for the getters which are not automatically computed
template <typename Type, typename PType>
Type CubicOrder<Type, PType> :: get_class_number(){
  if (class_number != 0){
    return class_number;
  } else{
    set_class_number();
    return class_number();
  }
}


template <typename Type, typename PType>
PType CubicOrder<Type, PType> :: get_regulator(){
  if (regulator != 0){
    return regulator;
  } else{
    set_regulator();
    return regulator();
  }
}


template <typename Type, typename PType>
std::vector<Type> CubicOrder<Type, PType> :: get_class_group(){
  if (cg_structure.size() != 0){
    return cg_structure;
  } else{
    set_class_group();
    return cg_structure();
  }
}

#endif
