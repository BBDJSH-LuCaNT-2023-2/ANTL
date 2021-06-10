#ifndef ANTL_COMPLEX_CUBIC_NUMBER_FIELD_CPP
#define ANTL_COMPLEX_CUBIC_NUMBER_FIELD_CPP


#include "../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"


/*
template <typename Type, typename PType>
RealCubicNumberField<Type, PType> :: RealCubicNumberField( polynomial<Type> const &poly)
  : CubicNumberFieldField<Type,PType>::CubicNumberField(poly) {

}
*/
template <typename Type, typename PType>
ComplexCubicNumberField<Type, PType> :: ComplexCubicNumberField( polynomial<Type> const &poly)
  : CubicNumberField<Type,PType>::CubicNumberField(poly) {

    if (this->discriminant > 0){
      std::cout << "Error, ComplexCubicNF input must be a polynomial with negative discriminant," << std::endl;
      throw std::exception();
    }

    set_roots();
    set_integral_basis();
}

template <typename Type, typename PType>
void ComplexCubicNumberField<Type, PType> :: set_roots(){
  PType roots[3];
  int root_type = cardano<PType>(this->definingPolynomial, roots);

  real_root = roots[0];
  complex_root = (roots[1], roots[2]);
};

template <typename Type, typename PType>
void ComplexCubicNumberField<Type, PType> :: set_integral_basis(){
  this->basis[0] = PType(1);                                     // 1
  mul(this->basis[1],(this->definingPolynomial)[3], real_root );  // a*delta

  this->basis[2] = this->basis[1];                              // a*delta
  add(this->basis[2],this->basis[2], (this->definingPolynomial)[2] ); // a*delta + b
  mul(this->basis[2], this->basis[2], real_root);
};











#endif // guard
