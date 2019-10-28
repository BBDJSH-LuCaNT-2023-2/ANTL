#ifndef ANTL_CUBIC_ORDER_H
#define ANTL_CUBIC_ORDER_H

/**
 * @file CubicOrder.hpp
 * @author Randy Yee
 * @brief order in a cubic number field or function field.
 */

#include <sstream>
#include <string>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include "../Arithmetic/QQ.hpp"
#include "../common.hpp"

#include "GeneralTemplateFunctions.hpp"
#include "CubicOrderReal.hpp"
#include "Multiplication/IdealMultiplicationStrategy.hpp"
#include "Multiplication/MultiplyStrategyWilliams.hpp"
#include "VoronoiMethods.hpp"
#include "VoronoiComplex.hpp"
#include "VoronoiReal.hpp"
using boost::math::tools::polynomial;
using NTL::abs;
// forward declaration
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class CubicIdeal;
template<typename Type, typename PType>
class IdealMultiplicationStrategy;
template<typename Type, typename PType>
class MultiplyStrategyWilliams;
template<typename Type, typename PType>
class VoronoiMethods;
template<typename Type, typename PType>
class VoronoiComplex;
template<typename Type, typename PType>
class VoronoiReal;

template <typename T, typename PT>
void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);
template <typename T, typename PT>
bool is_equal(CubicIdeal <T,PT> & A, CubicIdeal <T,PT> & B);

template <typename Type, typename PType>
class CubicOrder {

public:

static CubicOrder * make_order(polynomial<Type> const &poly){
  if (discriminant_bcf(poly) > 0 ){
    return new CubicOrderReal<Type, PType>(poly);
  }
  else{
    return new CubicOrder(poly);
  }
}
/** Constructor
* @param a polynomial to be used as the defining polynomial of the order.
* For NFs, an irreducible binary cubic form corresponds to a cubic order via Davenport-Heilbronn
*/
CubicOrder(polynomial<Type> const &poly, std::string s = "Williams");

/**
* @brief Constructor which requires a defining polynomial as well as a multiplication strategy
* @param A cubic polynomial, a pointer to a IdealMultiplicationStrategy object
*/
//CubicOrderNF<Type,PType> (polynomial<Type> const &poly, IdealMultiplicationStrategy<Type, PType> * mulstrat);



// Accessor methods
boost::math::tools::polynomial<Type> get_IBCF() const {return defining_IBCF;}
Type get_coeff(int i) const {
    return defining_IBCF[i];
}
std::shared_ptr<VoronoiMethods<Type, PType>> get_voronoi() const{
  return vmethods;
}

Type get_discriminant() const {return discriminant;}
PType get_rho1() const {return rho1;}
PType get_rho2() const {return rho2;}
PType get_root1() const {return root_list[0];}



// This is bad form, because in real these are the conjugates
// but in complex, this is the real and complex part respectively
inline PType get_root2() const {return root_list[1];}
inline PType get_root3() const {return root_list[2];}

inline long get_index() const {return index;}

ZZ get_class_number();
PType get_regulator();
std::vector<Type> get_class_group();

virtual CubicElement<Type, PType> * get_fund_unit(int i = 0){
  if (this->fundamentalUnits.size() == 0){
    this->compute_fundamental_unit();
  }
  return &fundamentalUnits[0];
};

inline bool is_real() const {
  return (discriminant > 0);
}

inline bool is_complex() const {
  return (discriminant < 0);
}
/**
* @brief Method for switching the order of the roots of the defining polynomial
* @param two integers which indicate the positions to swap
*/
virtual void roots_swap_position(int p1, int p2);

/**
* @brief Method for determining the splitting type of a rational prime
* @param A prime number p
*/
int splitting_type(int p);

/**
* @brief Method for setting the IdealMultiplicationStrategy
* @param A pointer to a concrete strategy (subclass of IdealMultiplicationStrategy)
*/
inline void set_mul_strategy(std::shared_ptr<IdealMultiplicationStrategy<Type,PType>> mstrat){
  this->m_method = std::static_pointer_cast< IdealMultiplicationStrategy<Type, PType>>(mstrat);
};



/** @brief returns the index of the order wrt the maximal order, the number f such that disc(O) = f^2 * Delta,
Delta = field discriminant.
* Currently no function has been implemented to compute the index. Use Llorente-Nart for NFs */
bool is_maximal(){if (index ==1) return true; else return false;}

/**
* @brief Method for checking if this CubicOrder is equal to a second one CO2
* @param A CubicOrder object)
*/
bool is_equal(const CubicOrder<Type, PType> &CO2) const;


/**
* @brief function to obtain the real value of a basis element, or conjugates
* @param[out]  newVal is the real value of the basis element
* @param[in] newVal is a float type, U,X,Y is the representation in terms of
*             {1, rho1, rho2}. In the real case can use a user-specified
*             conjugate of the defining root 0, 1, or 2.
*/
virtual void get_real_value(PType & newVal, Type &U, Type &X, Type &Y, Type &D, int conj = 0);





// **************************************************************************  /
// ******************* Friend classes and functions *************************  /

template <typename T, typename PT>
friend class IdealMultiplicationStrategy;
/** Can I delete this? I wanted mul is implemented in CubicOrderNF, but needed
* access to variables of CubicOrder. */
template <typename T, typename PT>
friend void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);

template <typename T, typename PT>
friend void mul(CubicIdeal<T,PT> & A, const CubicIdeal<T,PT> & B, const CubicIdeal<T,PT> & C);

template <typename T, typename PT>
friend bool is_equal(const CubicOrder<T, PT> &CO1, const CubicOrder<T, PT> &CO2);

protected:

// GlobalCubicField * my_cubic_field;
boost::math::tools::polynomial<Type> defining_IBCF;
//IdealMultiplicationStrategy<Type,PType> * mul_method;
std::shared_ptr<IdealMultiplicationStrategy<Type,PType>> m_method;
// Put pointers to concrete Multiply strategies here
std::shared_ptr<MultiplyStrategyWilliams<Type, PType>>  williams;
std::shared_ptr<VoronoiComplex<Type, PType>>  vcomplex;
std::shared_ptr<VoronoiReal<Type, PType>>  vreal;

////////////////////////////////////////////////////////////////

std::shared_ptr<VoronoiMethods<Type, PType>>  vmethods;

PType root_list[3];

Type discriminant;

// 2nd and 3rd elements of the integral basis
PType rho1, rho2;

long index; // The number f such that disc(O) = f^2 * Delta, Delta = field discriminant

Type class_number = Type(0);

std::vector<Type> cg_structure;

PType regulator = PType(0);

Type mul_table[3][3];

std::vector<CubicElement<Type, PType>> fundamentalUnits;

static PType order_temp;

CubicElement<Type, PType> temp_element = CubicElement<Type, PType>(this, Type(1), Type(0), Type(0), Type(1));

void set_roots();

// do I even need this?
void set_mul_table();

void set_integral_basis();

void set_class_number();

void set_class_group();

void set_regulator();

void compute_fundamental_unit();

private:

CubicOrder(); // private default constructor means you can't call this




}; // closes class scope


template<typename T, typename PT> PT CubicOrder<T,PT>::order_temp;
//template<typename T, typename PT> CubicElement<T, PT> CubicOrder<T,PT>::temp_element(T(1), T(0), T(0), T(1));
#include "../../../src/Cubic/CubicOrder.cpp"

#endif // guard
