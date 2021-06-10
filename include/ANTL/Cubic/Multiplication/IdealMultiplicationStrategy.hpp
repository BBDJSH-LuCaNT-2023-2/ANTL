#ifndef ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP
#define ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP

// strategy class for ideal multiplications
#include "../CubicOrder.hpp"


// forward declaration
template<typename Type, typename PType>
class CubicOrder;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class CubicIdeal;


template<typename Type, typename PType>
class IdealMultiplicationStrategy {

public:

  ~IdealMultiplicationStrategy(){

  };
    // this method should do some checking to ensure that the ideals are in the same order



//this will be a virtual function which is instantiated in subclasses
// different subclasses shall implement distinct methods
virtual void multiply(CubicIdeal<Type,PType> &A, const CubicIdeal<Type,PType> &B, const CubicIdeal<Type,PType> &C) = 0;


protected:

const CubicOrder<Type, PType> * my_order = NULL;
// variables to hold the product of two ideals as an integral 3 by 3 matrix along with a denominator.
Type product_basis[3][3];
Type multemp, denom;


// multiplies ideal basis elements internally
void column_mul(Type & out1, Type & out2, Type & out3, const Type & u1, const Type & x1, const Type & y1,
const Type & u2, const Type & x2, const Type & y2);
private:



}; //end class def

#include "../../../../src/Cubic/Multiplication/IdealMultiplicationStrategy.cpp"

#endif // include guard
