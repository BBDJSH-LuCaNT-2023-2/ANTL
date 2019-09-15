#ifndef ANTL_CUBIC_IDEAL_NF_CPP
#define ANTL_CUBIC_IDEAL_NF_CPP

#include "../../include/ANTL/Cubic/CubicIdealNF.hpp"


template <typename Type, typename PType>
CubicIdealNF<Type,PType> :: CubicIdealNF(const CubicOrderNF<Type,PType> * cnfo, const CubicElementNF<Type,PType> & A,
  const CubicElementNF<Type,PType> & B, const CubicElementNF<Type,PType> & C) : CubicIdeal<Type,PType>::CubicIdeal(cnfo),
  gen1(cnfo, Type(1),Type(0),Type(0),Type(1)),
  gen2(cnfo, Type(1),Type(0),Type(0),Type(1)),
  gen3(cnfo, Type(1),Type(0),Type(0),Type(1)) {
    gen1.assign(A);
    gen2.assign(B);
    gen3.assign(C);
    normalize();
  }


template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: normalize(){

  //set the new denominator
  mul(this->ci_temp, gen1.get_denom(), gen2.get_denom());
  mul(this->ci_temp, this->ci_temp, gen3.get_denom()); //ci_temp = d1*d2*d3
  mul(this->denom, this->denom, this->ci_temp);

  // Note we should be careful here since mul automatically normalizes after
  // completing the operation. By multiplying by d1*d2*d3, the denominators
  // of each generator should be 1. so there should be no problems normalizing
  // ideal
  mul(gen1, gen1, this->ci_temp);
  mul(gen2, gen2, this->ci_temp);
  mul(gen3, gen3, this->ci_temp);
  this->ci_temp = GCD(GCD(gen1.get_u(), gen2.get_u()),gen3.get_u() );

  if (this->ci_temp != Type(1)){
    this->ci_temp = GCD(GCD(GCD(GCD(GCD(GCD(GCD(
    this->ci_temp,
    gen1.get_x() ),
    gen2.get_x() ),
    gen3.get_x() ),
    gen1.get_y() ),
    gen2.get_y() ),
    gen3.get_y() ),
    this->denom  );

    div(gen1, gen1, this->ci_temp);
    div(gen2, gen2, this->ci_temp);
    div(gen3, gen3, this->ci_temp);
    div(this->denom, this->denom, this->ci_temp);

  }

}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_equivalent(const CubicIdeal<Type, PType> & B){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
Type CubicIdealNF<Type, PType> :: norm(){
  std::cout <<"undefined" << std::endl;
  return Type(-1);
}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_principal(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_prime(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_canonical(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: reduce(){
  std::cout <<"undefined" << std::endl;
}

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: become_canonical(){
  std::cout <<"undefined" << std::endl;
}

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: become_prepared(){
  std::cout <<"undefined" << std::endl;
}

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: voronoi(){
  std::cout <<"undefined" << std::endl;
}
#endif
