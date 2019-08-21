#ifndef ANTL_CUBIC_NUMBER_FIELD_CPP
#define ANTL_CUBIC_NUMBER_FIELD_CPP


#include <NTL/ZZX.h>
#include "../../include/ANTL/Cubic/CubicNumberField.hpp"

using NTL::ZZX;

/*
template <typename Type, typename PType>
CubicNumberField<Type, PType> *CubicNumberField<Type, PType>::create (polynomial<Type> const &poly) {
  Type disc = poly[2]*poly[2]*poly[1]*poly[1]
      + Type(18)*poly[3]*poly[2]*poly[1]*poly[0]
      - Type(4)*poly[3]*poly[1]*poly[1]*poly[1]
      - Type(4)*poly[2]*poly[2]*poly[2]*poly[0]
      - Type(27)*poly[3]*poly[3]*poly[0]*poly[0] ;

  return new RealCubicNumberField<long,long>;
}
*/

template <typename Type, typename PType>
CubicNumberField<Type, PType> :: CubicNumberField( polynomial<Type> const &poly)
  : GlobalCubicField<Type,PType>::GlobalCubicField(poly) {
    calc_discriminant();
}

template <typename Type, typename PType>
CubicNumberField<Type, PType> :: ~CubicNumberField (){

}






#endif // guard
