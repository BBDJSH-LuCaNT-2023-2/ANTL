#ifndef ANTL_COMPLEX_CUBIC_NUMBER_FIELD_H
#define ANTL_COMPLEX_CUBIC_NUMBER_FIELD_H

#include <NTL/ZZX.h>
#include "GlobalCubicField.hpp"
#include "CubicNumberField.hpp"



template <typename Type, typename PType>
class ComplexCubicNumberField : public CubicNumberField<Type, PType>{

public:


  // constructors
  ComplexCubicNumberField(polynomial<Type> const &poly);


  // setters
  void set_roots();

  void set_integral_basis();

  // getters
  PType get_real_root(){return real_root;}
  std::complex<PType>  get_complex_root(){return complex_root;}

private:

  PType real_root;
  std::complex<PType> complex_root;
}; //close class definition




#include "../../../src/Cubic/ComplexCubicNumberField.cpp"

#endif // guard
