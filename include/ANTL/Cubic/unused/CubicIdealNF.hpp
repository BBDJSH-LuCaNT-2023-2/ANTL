#ifndef ANTL_CUBIC_IDEAL_NF_H
#define ANTL_CUBIC_IDEAL_NF_H

/**
 * @file CubicIdealNF.hpp
 * @author Randy Yee
 * @remarks Class representing ideals of cubic number fields.
 */



#include "CubicOrderNF.hpp"
#include "CubicElement.hpp"
#include "CubicIdeal.hpp"
#include "../Arithmetic/QQ.hpp"

using namespace ANTL;

// class forward declarations
template<typename Type,typename PType>
class CubicOrderNF;
template<typename Type,typename PType>
class CubicElement;
template<typename Type,typename PType>
class CubicIdeal;
template<typename Type,typename PType>
class IdealMultiplicationStrategy;


template<typename Type,typename PType>
class CubicIdealNF : public CubicIdeal<Type, PType>{

template<typename T,typename PT>
friend class IdealMultiplicationStrategy;


public:


/************** Constructor(s) **********************/
CubicIdealNF(const CubicOrderNF<Type,PType> * cnfo, const CubicElement<Type,PType> & A,
  const CubicElement<Type,PType> & B, const CubicElement<Type,PType> & C);


/************** Accessors **********************/

inline const CubicElement<Type, PType> * get_gen1() const {return &gen1;}
inline const CubicElement<Type, PType> * get_gen2() const {return &gen2;}
inline const CubicElement<Type, PType> * get_gen3() const {return &gen3;}
inline const Type get_denom () const {return denom;}



void assign(const CubicElement<Type, PType> g3, const CubicElement<Type, PType> g2, const CubicElement<Type, PType> g1);
void assign(const Type U1, const Type X1, const Type Y1, const Type D1,
const Type U2, const Type X2, const Type Y2, const Type D2,
const Type U3, const Type X3, const Type Y3, const Type D3);
//


bool is_integral(){
  return IsOne(this->denom);
};

/**
* @brief this function checks whether two ideals are the same by comparing their
* lattices
*/
bool is_equal(const CubicIdeal<Type, PType> & B);

/**
* @brief function to determine whether the ideal is equivalent to B. Not yet implemented
*/
bool is_equivalent(const CubicIdeal<Type, PType> & B);


Type norm();
/**
* @brief function to determine whether the ideal is principal. Not yet implemented
*/
bool is_principal();

/**
* @brief function to determine whether the ideal is prime. Not yet implemented
*/
bool is_prime();

/**
* @brief function to determine whether the ideal is in canonical form. Not yet implemented
*/

bool is_canonical();

/**
* @brief function to obtain an equivalent reduced ideal. Not yet implemented
*/
void reduce();





/******************** friends *********************/
/**
* @brief the mul method
*
*/
//template <typename T, typename PT>
//friend void MultiplyStrategyWilliams<T, PT> :: multiply(CubicIdealNF<Type,PType>  &A, const CubicIdealNF<Type,PType>  &B, const CubicIdealNF<Type,PType> &C);


protected:

  // a Z basis for the ideal, should be of size 3
  // The representation will always be as 3 integral elements all over a common denominator
  static Type ci_temp, ci_temp2, ci_temp3;
  CubicElement<Type, PType> gen1;
  CubicElement<Type, PType> gen2;
  CubicElement<Type, PType> gen3;
  Type denom = Type(1);

  void normalize();



private:




}; // close class definition

template<typename T, typename PT> T CubicIdealNF<T,PT>::ci_temp;
template<typename T, typename PT> T CubicIdealNF<T,PT>::ci_temp2;
template<typename T, typename PT> T CubicIdealNF<T,PT>::ci_temp3;


#include "../../../src/Cubic/CubicIdealNF.cpp"
#endif
