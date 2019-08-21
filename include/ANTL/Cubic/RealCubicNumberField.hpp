#ifndef ANTL_REAL_CUBIC_NUMBER_FIELD_H
#define ANTL_REAL_CUBIC_NUMBER_FIELD_H

#include <NTL/ZZX.h>
#include "GlobalCubicField.hpp"
#include "CubicNumberField.hpp"
#include "GeneralTemplateFunctions.hpp"


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

  virtual void set_integral_basis();  // override
private:

  PType roots[3];


}; //close class definition

#include "../../../src/Cubic/RealCubicNumberField.cpp"

#endif // guard
