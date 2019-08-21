#ifndef ANTL_GLOBAL_CUBIC_FIELD_CPP
#define ANTL_GLOBAL_CUBIC_FIELD_CPP

#include "../../include/ANTL/Cubic/GlobalCubicField.hpp"
#include <boost/multiprecision/gmp.hpp>

using boost::multiprecision::mpf_float_100;


template <typename Type, typename PType>
GlobalCubicField<Type, PType>::GlobalCubicField(Type &A, Type &B, Type &C, Type &D) {
    definingPolynomial[3] = A;
    definingPolynomial[2] = B;
    definingPolynomial[1] = C;
    definingPolynomial[0] = D;
}

template <typename Type, typename PType>
GlobalCubicField<Type, PType>::GlobalCubicField(polynomial<Type> const &poly) {
  definingPolynomial = poly;

}

template <typename Type, typename PType>
GlobalCubicField<Type, PType>::GlobalCubicField(long A, long B, long C, long D) {
  definingPolynomial[3] = A;
  definingPolynomial[2] = B;
  definingPolynomial[1] = C;
  definingPolynomial[0] = D;
}

template <typename Type, typename PType>
GlobalCubicField<Type, PType>::GlobalCubicField() {
    definingPolynomial[3] = Type(1);
}


template <typename Type, typename PType>
GlobalCubicField<Type, PType> :: ~GlobalCubicField (){

}


#endif // guard
