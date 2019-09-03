#ifndef ANTL_REAL_CUBIC_NUMBER_FIELD_H
#define ANTL_REAL_CUBIC_NUMBER_FIELD_H

#include "CubicNumberField.hpp"


template <typename Type, typename PType>
class RealCubicNumberField : public CubicNumberField<Type, PType>{

public:

  RealCubicNumberField(polynomial<Type> const &poly);

  void roots_swap_position(int p1, int p2);

  inline PType get_root(int i){
    return roots[i];
  };

protected:
  void set_roots();

  /**
  * @brief Function which computes the canonical integral basis {1, a*rho1, a*rho^2 +b}
  * @pre The assumption is that the definingPolynomial has equation order 1,
  * otherwise we get the integral basis for some sub-order.
  */
  void set_integral_basis();  // override
private:

  PType roots[3];


}; //close class definition

#include "../../../src/Cubic/RealCubicNumberField.cpp"

#endif // guard
