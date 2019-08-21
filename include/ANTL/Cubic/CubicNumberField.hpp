#ifndef ANTL_CUBIC_NUMBER_FIELD_H
#define ANTL_CUBIC_NUMBER_FIELD_H


#include <NTL/ZZX.h>
#include "GlobalCubicField.hpp"
//#include "RealCubicNumberField.hpp"


template <typename Type, typename PType>
class CubicNumberField : public GlobalCubicField<Type, PType>{

public:

  // static factory method... is this the best way to do this?
  static CubicNumberField<Type,PType> *create (polynomial<Type> const &poly);

  // constructors and destructor
  CubicNumberField(polynomial<Type> const &poly);

  ~CubicNumberField ();

  bool is_complex(){ if (this->discriminant < Type(0)) {return true;} else {return false;}};
  bool is_real(){ if (this->discriminant > Type(0)) {return true;} else {return false;}};


protected:


  /**
  * uses the defining polynomial to compute the discriminant Delta = b^2c^2 + 18abcd - 4ac^3-rb^3d -27a^2d^2
  * Here, definingPolynomial[i] is the coefficient of the ith power of x (see boost documentation)
  */
  void calc_discriminant(){

    this->discriminant = this->definingPolynomial[2]*this->definingPolynomial[2]*this->definingPolynomial[1]*this->definingPolynomial[1]
    + Type(18)*this->definingPolynomial[3]*this->definingPolynomial[2]*this->definingPolynomial[1]*this->definingPolynomial[0]
    - Type(4)*this->definingPolynomial[3]*this->definingPolynomial[1]*this->definingPolynomial[1]*this->definingPolynomial[1]
    - Type(4)*this->definingPolynomial[2]*this->definingPolynomial[2]*this->definingPolynomial[2]*this->definingPolynomial[0]
    - Type(27)*this->definingPolynomial[3]*this->definingPolynomial[3]*this->definingPolynomial[0]*this->definingPolynomial[0] ;

  };  // override


  // Will use different data types and method for finding the roots
  // https://stackoverflow.com/questions/13328676/c-solving-cubic-equations
  virtual void set_roots() = 0;
}; //close class definition


#include "../../../src/Cubic/CubicNumberField.cpp"
#include "../../../src/Cubic/GlobalCubicField.cpp"

#endif // guard
