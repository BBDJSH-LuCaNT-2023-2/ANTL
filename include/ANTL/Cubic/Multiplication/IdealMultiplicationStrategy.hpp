#ifndef ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP
#define ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP

// strategy class for ideal multiplications
#include "../CubicIdeal.hpp"
#include "../CubicOrderNF.hpp"

// forward declaration
template<typename Type, typename PType>
class CubicOrderNF;


template<typename Type, typename PType>
class IdealMultiplicationStrategy {

public:


    // this method should do some checking to ensure that the ideals are in the same order



//this will be a virtual function which is instantiated in subclasses
// different subclasses shall implement distinct methods
//virtual void multiply(CubicIdeal<Type,PType> &C, const CubicIdeal<Type,PType> &A, const CubicIdeal<Type,PType> &B) = 0;



protected:

const CubicOrderNF<Type, PType> * my_order;
// variables to hold the product of two ideals as an integral 3 by 3 matrix along with a denominator.
Type product_basis[3][3];
Type denom;

private:



}; //end class def

#endif // include guard
