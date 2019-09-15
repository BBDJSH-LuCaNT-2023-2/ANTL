#ifndef ANTL_CUBIC_IDEAL_H
#define ANTL_CUBIC_IDEAL_H

/**
 * @file CubicIdeal.hpp
 * @author Randy Yee
 * @remarks Class representing a generic cubic ideal.
 */


#include "CubicOrder.hpp"
#include "CubicElement.hpp"
#include "CubicOrderNF.hpp"
#include "CubicElementNF.hpp"
#include "../Arithmetic/QQ.hpp"

using namespace ANTL;

//forward declarations
template<typename Type, typename PType>
class CubicOrder;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class MultiplyStrategyWilliams;


template<typename Type,typename PType>
class CubicIdeal{


public:

CubicIdeal(const CubicOrder<Type,PType> * cnfo){
  this->my_order = cnfo;
}

/************** Accessors **********************/
inline const CubicOrder<Type, PType> * get_order() const {return my_order;}



virtual bool is_equivalent(const CubicIdeal<Type, PType> & B);
virtual Type norm();
//

/******************** friends *********************/
//friend void mul<Type, PType > (CubicIdeal<Type, PType> & A, const CubicIdeal<Type, PType> & B, const CubicIdeal<Type, PType> & C);

protected:


// holds a reference to the order in which the ideal sits in
const CubicOrder<Type, PType> * my_order;


virtual void normalize();


private:




}; // close class definition
#include "../../../src/Cubic/CubicIdeal.cpp"
#endif
