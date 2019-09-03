#ifndef ANTL_CUBIC_ELEMENT_CPP
#define ANTL_CUBIC_ELEMENT_CPP

#include "../../include/ANTL/Cubic/CubicElement.hpp"


template<typename Type, typename PType>
Type CubicElement<Type, PType>::temp;
template<typename Type, typename PType>
Type CubicElement<Type, PType>::newU;
template<typename Type, typename PType>
Type CubicElement<Type, PType>::newX;
template<typename Type, typename PType>
Type CubicElement<Type, PType>::newY;

//long CubicElement<long, float>::temp = 0;


template <typename Type, typename PType>
void CubicElement<Type, PType> :: assign(const CubicElement<Type,PType> & C){
  this->u = C.get_u();
  this->x = C.get_x();
  this->y = C.get_y();
  this->denom = C.get_denom();

  this->normalize();

}

template <typename Type, typename PType>
void CubicElement<Type, PType> :: assign(const Type _coeff[3], const Type & D) {
  this->u = _coeff[0];
  this->x = _coeff[1];
  this->y = _coeff[2];
  this->denom = D;

  this->normalize();
};

template <typename Type, typename PType>
void CubicElement<Type, PType> :: assign(const Type & U,const Type & X,const Type & Y, const Type & D) {
  this->u = U;
  this->x = X;
  this->y = Y;
  this->denom = D;

  this->normalize();

};

template <typename Type, typename PType>
bool CubicElement<Type, PType> :: is_equal(const CubicElement <Type, PType> & B){
  return ( ( (this->my_order)->is_equal(*B.get_order()) ) && (this->u == B.u) && (this->x == B.x) && (this->y == B.y)&& (this->denom == B.denom));
}
#endif
