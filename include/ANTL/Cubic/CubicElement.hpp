#ifndef ANTL_CUBIC_ELEMENT_H
#define ANTL_CUBIC_ELEMENT_H

#include "CubicOrder.hpp"

using namespace NTL;
//#include "../Arithmetic/QQ.hpp"



template<typename Type, typename PType>
class CubicElement{

public:

/***************** member variables **********************/

/***************** member functions **********************/

CubicElement(const CubicOrder<Type,PType> * cnfo, const Type _coefficients[3], const Type & _denom){
  u = _coefficients[0];
  x = _coefficients[1];
  y = _coefficients[2];
  denom = _denom;
  my_order = cnfo;
}



/*********** Accessor functions **************/
inline const CubicOrder<Type, PType> * get_order() const {return my_order;}

inline Type get_denom () const {return denom;}
inline Type get_u () const {return u;}
inline Type get_x () const {return x;}
inline Type get_y () const {return y;}


/******calculator functions****/
//these should calculate and set into the provided variable the specified value

//virtual void norm(ANTL::QQ<Type> & newVal) = 0;
//virtual void trace(ANTL::QQ<Type> & newVal) = 0;

/** @brief Should set this CubicElement to C */
void assign(const CubicElement<Type,PType> & C) {
  this->u = C.u;
  this->x = C.x;
  this->y = C.y;
  this->denom = C.denom;

};

//this should specify whether val is zero
virtual bool is_zero() const = 0;

//these procedural operations should take B and C and place the result into A



// The functions below turned into friend functions, no need for member operations?
//void subtract (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const CubicElement <Type,PType>& C);
//void divide (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const CubicElement <Type,PType> & C);
//void power (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const NTL::ZZ & p);




//template <typename T, typename PT>
//friend void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);


protected:


/***************** member variables **********************/


const CubicOrder<Type, PType> * my_order; /** A reference to the order which the element belongs to. */

Type u,x,y;           /** Coefficients in terms of an integral basis of my_order. It is understood that alpha = (u + x*rho1 + y*rho2)/denom */
Type denom;           /** Common denominator of coefficients */

/** Temporary variable(s) for arithmetic operations */
static Type newU, newX, newY,temp;
/***************** member functions **********************/


virtual void reduce()=0; // This function should reduce the coefficient vector and denominator into lowest terms.

private:


}; // close class definition




// static forward declaration

#include "../../../src/Cubic/CubicElement.cpp"

#endif
