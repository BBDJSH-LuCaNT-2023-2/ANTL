#ifndef ANTL_CUBIC_ELEMENT_H
#define ANTL_CUBIC_ELEMENT_H

#include "CubicOrder.hpp"
#include <boost/multiprecision/gmp.hpp>
#include "../Arithmetic/QQ.hpp"
using namespace ANTL;
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
* This is kind of useless since we haven't assigned a cubic order
* Any element initialized with the default constructor must be assigned an order
*/
CubicElement(){
  NTL::set(this->denom);
  NTL::set(this->u);
  NTL::clear(this->x);
  NTL::clear(this->y);
};

/**
* @brief Constructor which takes in a cubic order and creates the one element
*/
CubicElement(const CubicOrder<Type,PType> * cnfo){
  my_order = cnfo;
  NTL::set(this->denom);
  NTL::set(this->u);
  NTL::clear(this->x);
  NTL::clear(this->y);
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
  NTL::set(this->denom);
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


inline bool is_zero() const {
  return (this->u ==Type(0)) && (this->x ==Type(0)) && (this->y == Type(0));
}
/**
* @brief Computes the inverse element of this CubicElement and puts it in val
* @ param[in] val a CubicElement
* @pre val must be in the same order as this CubicElement
*/

void inverse(CubicElement<Type, PType> & val) const;
//this should specify whether val is zero




/**
* @brief Computes the absolute norm of this CubicElement
* and puts it in newVal
*/
void norm(ANTL::QQ<Type> & newVal);
/**
* @brief Computes the absolute trace of this CubicElement
* and puts it in newVal
*/
void trace(ANTL::QQ<Type> & newVal);


/**
* @brief Return the real value of the CubicElement
* Obtained by subbing in the main root of the associated CubicOrder
*/
void get_real_value(PType & newVal) ;


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
static Type newU, newX, newY,temp, temp2;
static PType precise_temp;
/***************** member functions **********************/


void normalize(); // This function should reduce the coefficient vector and denominator into lowest terms.


private:



  /************************** friend functions ********************************/
  /************************ forward declarations ******************************/

  /**
  * @brief Sets A to equal B
  * @param[out] A is the result of assignment
  * @param[in] B the value to be assigned
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void assign (CubicElement <T,PT> & A, const CubicElement <T,PT> & B);

  /**
  * @brief Compute product of B and C
  * @param[out] A = B * C
  * @param[in] B first term
  * @param[in] C second term
  * @pre A, B, C should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);

  /**
  * @brief Compute product of B and alpha (an integer constant)
  * @param[out] A = B * alpha
  * @param[in] B first term
  * @param[in] alpha second term (constant)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const T & alpha);

  /**
  * @brief Compute product of B and alpha (rational number)
  * @param[out] A = B * alpha
  * @param[in] B first term
  * @param[in] alpha second term (rational number)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const QQ<T> & alpha);


  /**
  * @brief Compute sum of B and C
  * @param[out] A = B + C
  * @param[in] B first term
  * @param[in] C second term
  * @pre A, B, C should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void add (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);

  /**
  * @brief Compute sum of B and alpha (an integer constant)
  * @param[out] A = B + alpha
  * @param[in] B first term
  * @param[in] alpha second term (constant)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void add (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const T & alpha);

  /**
  * @brief Compute sum of B and alpha (rational number)
  * @param[out] A = B + alpha
  * @param[in] B first term
  * @param[in] alpha second term (rational number)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void add (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const QQ<T> & alpha);

  /**
  * @brief Compute difference of B and C
  * @param[out] A = B - C
  * @param[in] B first term
  * @param[in] C second term
  * @pre A, B, C should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void sub (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT>& C);

  /**
  * @brief Compute difference of B and alpha (an integer  constant)
  * @param[out] A = B - alpha
  * @param[in] B first term
  * @param[in] alpha second term (constant)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void sub (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const T & alpha);

  /**
  * @brief Compute difference of B and alpha (rational number)
  * @param[out] A = B - alpha
  * @param[in] B first term
  * @param[in] alpha second term (rational number)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void sub (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const QQ<T> & alpha);


  /**
  * @brief Compute quotient B / C where B,C are both cubic elements. Not yet implemented
  * @param[out] A = B / Ca
  * @param[in] B first term
  * @param[in] C second term
  * @pre A, B, C should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void div (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);
  /**
  * @brief Compute quotient of B by alpha (an integer constant)
  * @param[out] A = B / alpha
  * @param[in] B first term
  * @param[in] alpha second term (constant)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void div (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const T & alpha);

  /**
  * @brief Compute quotient of B and alpha (a rational number)
  * @param[out] A = B / alpha
  * @param[in] B first term
  * @param[in] alpha second term (rational number)
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void div (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const QQ<T> & alpha);
  /* We can't divide in the ring of integers, but I want to be able to invert field elements?
  template <typename T, typename PT>
  friend void divide (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT>& C);
  */



  /**
  * @brief NOT YET IMPLEMENTED Compute square of B.
  * @param[out] A = B*B
  * @param[in] B first term
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void sqr (CubicElement <T,PT> & A, const CubicElement <T,PT> & B);

  /**
  * @brief NOT YET IMPLEMENTED Compute cube of B.
  * @param[out] A = B*B*B
  * @param[in] B first term
  * @pre A, B should belong to the same CubicOrder
  */
  template <typename T, typename PT>
  friend void cube (CubicElement <T,PT> & A, const CubicElement <T,PT> & B);


}; // close class definition




// static forward declaration
#include "../../../src/Cubic/CubicElement.cpp"

#endif
