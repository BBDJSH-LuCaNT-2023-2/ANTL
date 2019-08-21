#ifndef ANTL_REAL_CUBIC_NUMBER_FIELD_CPP
#define ANTL_REAL_CUBIC_NUMBER_FIELD_CPP


#include <NTL/ZZX.h>
#include "../../include/ANTL/Cubic/RealCubicNumberField.hpp"


/*
template <typename Type, typename PType>
RealCubicNumberField<Type, PType> :: RealCubicNumberField( polynomial<Type> const &poly)
  : CubicNumberFieldField<Type,PType>::CubicNumberField(poly) {

}
*/
template <typename Type, typename PType>
RealCubicNumberField<Type, PType> :: RealCubicNumberField( polynomial<Type> const &poly)
  : CubicNumberField<Type,PType>::CubicNumberField(poly) {

    if (this->discriminant < 0){
      std::cout << "error, input must be a polynomial with positive discriminant " << std::endl;
      throw std::exception();
    }
    set_roots();

}

template <typename Type, typename PType>
void RealCubicNumberField<Type, PType> ::  roots_swap_position(int p1, int p2){
  std::swap(this->roots[p1], this->roots[p2]);

};

template <typename Type, typename PType>
void RealCubicNumberField<Type, PType> ::  set_roots(){
  PType a1 = this->definingPolynomial[2]/this->definingPolynomial[3];
  PType b1 = this->definingPolynomial[1]/this->definingPolynomial[3];
  PType c1 = this->definingPolynomial[0]/this->definingPolynomial[3];

  int root_type = SolveP3<PType>(this->roots, a1,b1,c1);
  if (abs(roots[0])-1.0 <= 0){
    roots_swap_position(this->roots[0], this->roots[1]);
  }
};

template <typename Type, typename PType>
void RealCubicNumberField<Type, PType> :: set_integral_basis(){
   this->basis[0] = PType(1);                                     // 1
   mul(this->basis[1],(this->definingPolynomial)[3], roots[0] );  // a*delta

   this->basis[2] = this->basis[1];                              // a*delta
   add(this->basis[2],this->basis[2], (this->definingPolynomial)[2] ); // a*delta + b
   mul(this->basis[2], this->basis[2], roots[0]);
 };



#endif // guard
