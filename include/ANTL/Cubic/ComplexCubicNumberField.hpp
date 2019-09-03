#ifndef ANTL_COMPLEX_CUBIC_NUMBER_FIELD_H
#define ANTL_COMPLEX_CUBIC_NUMBER_FIELD_H


#include "CubicNumberField.hpp"



template <typename Type, typename PType>
class ComplexCubicNumberField : public CubicNumberField<Type, PType>{

public:


  // constructors
  ComplexCubicNumberField(polynomial<Type> const &poly);

  // getters
  PType get_real_root(){return real_root;}
  std::complex<PType>  get_complex_root(){return complex_root;}

private:


// ********************* Data members ***************************//
  PType real_root;
  std::complex<PType> complex_root;



  // setters

  void set_roots();
  void set_integral_basis();

}; //close class definition




#include "../../../src/Cubic/ComplexCubicNumberField.cpp"

#endif // guard
