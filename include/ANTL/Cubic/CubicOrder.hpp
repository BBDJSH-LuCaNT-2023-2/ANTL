#ifndef ANTL_CUBIC_ORDER_H
#define ANTL_CUBIC_ORDER_H

/**
 * @file CubicOrder.hpp
 * @author Randy Yee
 * @brief order in a cubic number field or function field. Template class which
 * should specialize to CubicOrderNF, and later CubicOrderFF
 */

#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include "../Arithmetic/QQ.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "../common.hpp"

#include "GeneralTemplateFunctions.hpp"
#include "Multiplication/IdealMultiplicationStrategy.hpp"
using boost::math::tools::polynomial;

// forward declaration
template<typename Type, typename PType>
class CubicElementNF;

template <typename Type, typename PType>
class CubicOrder {

public:

/** Constructor
* @param a polynomial to be used as the defining polynomial of the order.
* For NFs, an irreducible binary cubic form corresponds to a cubic order via Davenport-Heilbronn
*/
CubicOrder(polynomial<Type> const &poly);




// Accessor methods
boost::math::tools::polynomial<Type> get_IBCF() const {return defining_IBCF;}
Type get_discriminant() const {return discriminant;}
Type get_rho1() const {return rho1;}
Type get_rho2() const {return rho2;}
long get_index() const {return index;}

Type get_class_number();
PType get_regulator();
std::vector<Type> get_class_group();


/** @brief returns the index of the order wrt the maximal order, the number f such that disc(O) = f^2 * Delta,
Delta = field discriminant.
* Currently no function has been implemented to compute the index. Use Llorente-Nart for NFs */
bool is_maximal(){if (index ==1) return true; else return false;}

/** Virtual function to test equality with another CubicOrder. */
virtual bool is_equal(const CubicOrder<Type, PType> &CO2) const = 0;



/***************** Friend functions ************************/

/** Can I delete this? I wanted mul is implemented in CubicOrderNF, but needed
* access to variables of CubicOrder. */
template <typename T, typename PT>
friend void mul (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT> & C);



protected:

// GlobalCubicField * my_cubic_field;
boost::math::tools::polynomial<Type> defining_IBCF;

IdealMultiplicationStrategy<Type, PType> * mul_method;

PType root_list[3];

Type discriminant;

// 2nd and 3rd elements of the integral basis
Type rho1, rho2;

long index; // The number f such that disc(O) = f^2 * Delta, Delta = field discriminant

Type class_number = Type(0);

std::vector<Type> cg_structure;

PType regulator = PType(0);


virtual void set_integral_basis() = 0;

virtual void set_class_number() = 0;

virtual void set_class_group() = 0;

virtual void set_regulator() = 0;





private:

CubicOrder(); // private default constructor means you can't call this




}; // closes class scope


#include "../../../src/Cubic/CubicOrder.cpp"

#endif // guard
