#ifndef ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP
#define ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP

// strategy class for ideal multiplications
#include "../CubicOrderNF.hpp"
//#include "../CubicElement.hpp"
//#include "../CubicIdeal.hpp"
//#include "../CubicIdealNF.hpp"

#include "../../XGCD/xgcd_plain.hpp"

// forward declaration
template<typename Type, typename PType>
class CubicOrderNF;
//template<typename Type, typename PType>
//class CubicElement;
//template<typename Type, typename PType>
//class CubicElementNF;
template<typename Type, typename PType>
class CubicIdealNF;


template<typename Type, typename PType>
class IdealMultiplicationStrategy {

public:

  ~IdealMultiplicationStrategy(){
    delete my_order;
  };
    // this method should do some checking to ensure that the ideals are in the same order



//this will be a virtual function which is instantiated in subclasses
// different subclasses shall implement distinct methods
virtual void multiply(CubicIdealNF<Type,PType> &C, const CubicIdealNF<Type,PType> &A, const CubicIdealNF<Type,PType> &B) = 0;


int testmember = 123;
protected:

const CubicOrderNF<Type, PType> * my_order;
// variables to hold the product of two ideals as an integral 3 by 3 matrix along with a denominator.
Type product_basis[3][3];
Type denom;

private:



}; //end class def

#endif // include guard
