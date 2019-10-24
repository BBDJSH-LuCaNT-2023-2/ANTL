#ifndef ANTL_CUBIC_NUMBER_FIELD_CPP
#define ANTL_CUBIC_NUMBER_FIELD_CPP


#include "../../include/ANTL/Cubic/CubicNumberField.hpp"



template <typename Type, typename PType>
CubicNumberField<Type, PType> :: CubicNumberField( polynomial<Type> const &poly)
  : GlobalCubicField<Type,PType>::GlobalCubicField(poly) {
    calc_discriminant();
}

template <typename Type, typename PType>
CubicNumberField<Type, PType> :: ~CubicNumberField (){

}






#endif // guard
