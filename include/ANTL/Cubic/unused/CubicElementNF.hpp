#ifndef ANTL_CUBIC_ELEMENTNF_H
#define ANTL_CUBIC_ELEMENTNF_H


#include "CubicElement.hpp"
#include "CubicOrderNF.hpp"

//forward declaration
//template<typename Type, typename PType>
//  class CubicOrder;
#include "../Arithmetic/QQ.hpp"
using namespace ANTL;

template<typename Type, typename PType>
class CubicOrderNF;


template <typename Type, typename PType>
class CubicElementNF : public CubicElement<Type, PType>{


public:

  /** Default constructor: I wanted this to be private,
  * but IdealMultiplicationStrategy needs it */
  CubicElementNF();


/**
* @brief Constructor which takes in all relevant values (order, coefficients, denominator)
* @param cnfo pointer to the CubicOrderNF which the element belongs
* @param _coefficients Values of the representation of the elemnent in terms of an integral basis
* @param denominator  Common denominator of the coefficients
*/
CubicElementNF(const CubicOrderNF<Type, PType> * cnfo, const Type _coefficients[3], const Type & _denom);
CubicElementNF(const CubicOrderNF<Type, PType> * cnfo, const Type U, const Type X, const Type Y, const Type & _denom);
//void norm(ANTL::QQ<Type> & newVal);
void trace(ANTL::QQ<Type> & newVal);
void inverse(CubicElementNF<Type,PType> & newVal){};

/**
* @brief Checks if this CubicElement is zero
*/
bool is_zero() const;




//void power (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const NTL::ZZ & p);

protected:

/**
* @brief Should put the coefficients and denominator of this CubicElementNF into
* lowest terms. Note that this is called normalize() in the quadratic library, consider changing it.
*/
void normalize();



private:
static Type temp2;






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
friend void mul (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT> & C);

/**
* @brief Compute product of B and alpha (an integer constant)
* @param[out] A = B * alpha
* @param[in] B first term
* @param[in] alpha second term (constant)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void mul (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const T & alpha);

/**
* @brief Compute product of B and alpha (rational number)
* @param[out] A = B * alpha
* @param[in] B first term
* @param[in] alpha second term (rational number)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void mul (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const QQ<T> & alpha);


/**
* @brief Compute sum of B and C
* @param[out] A = B + C
* @param[in] B first term
* @param[in] C second term
* @pre A, B, C should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void add (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT> & C);

/**
* @brief Compute sum of B and alpha (an integer constant)
* @param[out] A = B + alpha
* @param[in] B first term
* @param[in] alpha second term (constant)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void add (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const T & alpha);

/**
* @brief Compute sum of B and alpha (rational number)
* @param[out] A = B + alpha
* @param[in] B first term
* @param[in] alpha second term (rational number)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void add (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const QQ<T> & alpha);

/**
* @brief Compute difference of B and C
* @param[out] A = B - C
* @param[in] B first term
* @param[in] C second term
* @pre A, B, C should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void sub (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT>& C);

/**
* @brief Compute difference of B and alpha (an integer  constant)
* @param[out] A = B - alpha
* @param[in] B first term
* @param[in] alpha second term (constant)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void sub (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const T & alpha);

/**
* @brief Compute difference of B and alpha (rational number)
* @param[out] A = B - alpha
* @param[in] B first term
* @param[in] alpha second term (rational number)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void sub (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const QQ<T> & alpha);


/**
* @brief Compute quotient B / C where B,C are both cubic elements. Not yet implemented
* @param[out] A = B / Ca
* @param[in] B first term
* @param[in] C second term
* @pre A, B, C should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void div (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT> & C);
/**
* @brief Compute quotient of B by alpha (an integer constant)
* @param[out] A = B / alpha
* @param[in] B first term
* @param[in] alpha second term (constant)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void div (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const T & alpha);

/**
* @brief Compute quotient of B and alpha (a rational number)
* @param[out] A = B / alpha
* @param[in] B first term
* @param[in] alpha second term (rational number)
* @pre A, B should belong to the same CubicOrder
*/
template <typename T, typename PT>
friend void div (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const QQ<T> & alpha);
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


};// close class def




//template class CubicElementNF<long, boost::multiprecision::mpf_float_100>;
#include "../../../src/Cubic/CubicElementNF.cpp"

#endif
