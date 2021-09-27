#ifndef ANTL_FUND_UNIT_STRATEGY_HPP
#define ANTL_FUND_UNIT_STRATEGY_HPP

// strategy class for ideal multiplications
#include "../CubicOrder.hpp"
#include "../CubicIdeal.hpp"

// forward declaration
template<typename Type, typename PType>
class CubicOrder;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class CubicIdeal;


template<typename Type, typename PType>
class FundUnitStrategy {

public:

  ~FundUnitStrategy(){

  };
    // this method should do some checking to ensure that the ideals are in the same order



//this will be a virtual function which is instantiated in subclasses
// different subclasses shall implement distinct methods
virtual void compute( std::vector<CubicElement<Type, PType>> & unitvec, CubicOrder<Type, PType> * ord, bool realorder) {};


protected:

//const CubicOrder<Type, PType> * my_order = NULL;
// variables to hold the product of two ideals as an integral 3 by 3 matrix along with a denominator.

// multiplies ideal basis elements internally


private:



}; //end class def


#endif // include guard
