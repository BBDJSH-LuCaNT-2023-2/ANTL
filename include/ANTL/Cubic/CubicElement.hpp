#ifndef ANTL_CUBIC_ELEMENT_H
#define ANTL_CUBIC_ELEMENT_H

#include "CubicOrder.hpp"

#include <boost/multiprecision/gmp.hpp>
using namespace NTL;
//#include "../Arithmetic/QQ.hpp"


template<typename Type, typename PType>
class CubicOrder;

template<typename Type, typename PType>
class CubicElement{

public:

/***************** member variables **********************/

/***************** member functions **********************/

/**
* @brief Default Constructor
*/
CubicElement(){

};

CubicElement(const CubicOrder<Type,PType> * cnfo, const Type _coefficients[3], const Type & _denom){

  //throw exception if denom is zero
  if (_denom == Type(0)){
    std::cout << "error, 3rd argument cannot be 0 " << std::endl;
    throw std::exception();
  }

  // assign values to basis coefficients, denom, and order
  u = _coefficients[0];
  x = _coefficients[1];
  y = _coefficients[2];
  denom = _denom;
  my_order = cnfo;
}

CubicElement(const CubicOrder<Type,PType> * cnfo, const Type U, const Type X, const Type Y, const Type & _denom){

  //throw exception if denom is zero
  if (_denom == Type(0)){
    std::cout << "error, 3rd argument cannot be 0 " << std::endl;
    throw std::exception();
  }

  // assign values to basis coefficients, denom, and order
  u = U;
  x = X;
  y = Y;
  denom = _denom;
  my_order = cnfo;
}

/*********** Accessor functions **************/
inline const CubicOrder<Type, PType> * get_order() const {return my_order;}

inline Type get_denom () const {return denom;}
inline Type get_u () const {return u;}
inline Type get_x () const {return x;}
inline Type get_y () const {return y;}

inline void set_order(const CubicOrder<Type, PType> * ord){
  my_order = ord;
}
/******calculator functions****/
//these should calculate and set into the provided variable the specified value

//virtual void norm(ANTL::QQ<Type> & newVal) = 0;
//virtual void trace(ANTL::QQ<Type> & newVal) = 0;

/** @brief Should set this CubicElement to C */

/**
* @brief Sets this CubicElementNF equal to C
* Note that there is a friend version too
* @pre this CubicElementNF and C must be in the same CubicOrder
*/
void assign(const CubicElement<Type,PType> & C);


void assign(const Type _coeff[3], const Type & D);

void assign(const Type & U,const Type & X,const Type & Y, const Type & D);

/**
 * @brief Sets this CubicElement equal to n
 * @param[in] n value to give the CubicElementNF
 */
void assign(const Type & n){
  this->u = n;
  ::clear(this->x);
  ::clear(this->y);
  ::set(this->denom);
}

/**
 * @brief Sets this CubicElement equal to r
 * @param[in] r value to give the CubicElementNF
 */
void assign(const QQ<Type> & r){
  this->u = r.getNumerator();
  ::clear(this->x);
  ::clear(this->y);
  this->denom = r.getDenominator();
}

/**
* @brief Checks if this CubicElement is equal to B
* @ param[in] B a CubicElementNF
*/
bool is_equal(const CubicElement <Type, PType> & B);


//this should specify whether val is zero
virtual bool is_zero() const = 0;






//these procedural operations should take B and C and place the result into A



// The functions below turned into friend functions, no need for member operations?
//void subtract (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const CubicElement <Type,PType>& C);
//void divide (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const CubicElement <Type,PType> & C);
//void power (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const NTL::ZZ & p);






/**
* @brief Sets A to equal B
* @param[out] A is the result of assignment
* @param[in] B the value to be assigned
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void assign (CubicElement <T,PT> & A, const CubicElement <T,PT> & B);


/*****************************************************/
/*****************************************************/
protected:


/***************** member variables ******************/


const CubicOrder<Type, PType> * my_order;   // A reference to the order which the element belongs to.

Type u,x,y;           // Coefficients in terms of an integral basis of my_order. It is understood that alpha = (u + x*rho1 + y*rho2)/denom */
Type denom;           // Common denominator of coefficients

/** Temporary variable(s) for arithmetic operations */
static Type newU, newX, newY,temp;
/***************** member functions **********************/


virtual void normalize()=0; // This function should reduce the coefficient vector and denominator into lowest terms.

private:


}; // close class definition




// static forward declaration
#include "../../../src/Cubic/CubicElement.cpp"

#endif
