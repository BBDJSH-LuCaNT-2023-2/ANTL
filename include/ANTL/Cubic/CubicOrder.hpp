#ifndef ANTL_CUBIC_ORDER_H
#define ANTL_CUBIC_ORDER_H

//#define DEBUG 1
//#define DEBUGVORONOI 1
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
#include  "../Exponentiation/Exponentiation.hpp"
#include "../common.hpp"

#include "GeneralTemplateFunctions.hpp"
#include "CubicOrderReal.hpp"
#include "CubicOrderComplex.hpp"
#include "Multiplication/IdealMultiplicationStrategy.hpp"
#include "Multiplication/MultiplyStrategyWilliams.hpp"
#include "FundamentalUnits/FundUnitStrategy.hpp"
#include "FundamentalUnits/BasicVoronoi.hpp"
#include "FundamentalUnits/BSGSVoronoi.hpp"
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
template<typename Type, typename PType>
class FundUnitStrategy;
template<typename Type, typename PType>
class BasicVoronoi;
template<typename Type, typename PType>
class BSGSVoronoi;

template <typename T, typename PT>
void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);
template <typename T, typename PT>
bool is_equal(CubicIdeal <T,PT> & A, CubicIdeal <T,PT> & B);

template <typename Type, typename PType>
class CubicOrder {

public:

static CubicOrder * make_order(polynomial<Type> const &poly){
  if (discriminant_bcf(poly) >= 0 ){
    return new CubicOrderReal<Type, PType>(poly);
  }
  else {
    return new CubicOrderComplex<Type, PType>(poly);
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

/**
* @brief Returns a pointer to the associated VoronoiMethods object
*/
std::shared_ptr<VoronoiMethods<Type, PType>> get_voronoi() const{
  return vmethods;
}

Type get_discriminant() const {return discriminant;}
PType get_rho1() const {return rho1;}
PType get_rho2() const {return rho2;}
PType get_root1() const {return root_list[0];}



// This is bad form, because in real these are the conjugates
// but in complex, this is the real and complex part respectively
PType get_root2() const {return root_list[1];}
PType get_root3() const {return root_list[2];}

long get_index() const {return index;}

ZZ get_class_number();
PType get_regulator();
std::vector<Type> get_class_group();





/**
* @brief Method returns a bool indicating if the order is real (positive discriminant)
*/
bool is_real() const {return (discriminant > 0);}
/**
* @brief Method returns a bool indicating if the order is complex (negative discriminant)
*/
bool is_complex() const {return (discriminant < 0);}
/**
* @brief Method for switching the order of the roots of the defining polynomial
* Note that in the real case, this function is specialized to re-compute the conjugates of all basis elements
* @param two integers which indicate the positions to swap
*/
virtual void roots_swap_position(int p1, int p2);

/**
* @brief Method for determining the splitting type of a rational prime
* @param A prime number p
*/
int splitting_type(Type p);

/**
* @brief Method for setting the IdealMultiplicationStrategy
* @param A pointer to a concrete strategy (subclass of IdealMultiplicationStrategy)
*/
void set_mul_strategy(std::shared_ptr<IdealMultiplicationStrategy<Type,PType>> mstrat){
  this->m_method = std::static_pointer_cast< IdealMultiplicationStrategy<Type, PType>>(mstrat);
};
void set_unit_strategy(std::shared_ptr<FundUnitStrategy<Type,PType>> ustrat){
  this->unit_strat = std::static_pointer_cast< FundUnitStrategy<Type, PType>>(ustrat);
};

/**
* @brief Method for setting the FundUnitStrategy
* @param A pointer to a concrete strategy (subclass of FundUnitStrategy)
*/
void set_unit_strategy(std::string s){
  if (s.compare("Voronoi") == 0){
      if (!(this->voronoi_basic)){
          this->voronoi_basic.reset();
          this->voronoi_basic = std::make_shared<BasicVoronoi<Type, PType>>();
      }
      this->unit_strat = std::static_pointer_cast< FundUnitStrategy<Type, PType>>(this->voronoi_basic);
    }
  else if ( (s.compare("BSGS") == 0) && (this->is_complex()) ){
    if (!(this->voronoi_bsgs)){
        this->voronoi_bsgs.reset();
        this->voronoi_bsgs = std::make_shared<BSGSVoronoi<Type, PType>>();
    }
    this->unit_strat = std::static_pointer_cast< FundUnitStrategy<Type, PType>>(this->voronoi_bsgs);
  }
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
* Given a generating polynomial f, obtain its standard form, x^3 + Ex + G,
* where E,G are integers. See pg 22 CFG for formulas
*/
void standard_form(Type & E, Type& G);


/**
* @brief function to obtain the real value of a basis element, or conjugates
* @param[out]  newVal is the real value of the basis element
* @param[in] newVal is a float type, U,X,Y is the representation in terms of
*             {1, rho1, rho2}. In the real case can use a user-specified
*             conjugate of the defining root 0, 1, or 2.
*/
virtual void get_real_value(PType & newVal, const Type &U, const Type &X, const Type &Y, const Type &D, int conj = 0)=0;

/**
* @brief Method returns the fundamental unit.
* If the computation has not been performed, it computes and stores it first.
*/
virtual CubicElement<Type, PType> * get_fundamental_unit(int i = 0) = 0;

/**
* @brief given an ideal1, and a log vector vec1 (length r), return a reduced ideal J and a minima
* whose logarithm vector is close to vec1, and such that J = 1/minimum * ideal1
*/
virtual void close_minimum(CubicIdeal<Type, PType> & reduced_ideal1, \
  CubicElement<Type, PType> & minimum, CubicIdeal<Type, PType> & ideal1, std::vector<PType> & vec1) = 0;
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

std::shared_ptr<IdealMultiplicationStrategy<Type,PType>> m_method;
// Put pointers to concrete Multiply strategies here
std::shared_ptr<MultiplyStrategyWilliams<Type, PType>>  williams;
std::shared_ptr<VoronoiComplex<Type, PType>>  vcomplex;
std::shared_ptr<VoronoiReal<Type, PType>>  vreal;


std::shared_ptr<FundUnitStrategy<Type,PType>> unit_strat;
// Put pointers to concrete unit strategies here
std::shared_ptr<BasicVoronoi<Type, PType>>  voronoi_basic;
std::shared_ptr<BSGSVoronoi<Type, PType>>  voronoi_bsgs;
////////////////////////////////////////////////////////////////

std::shared_ptr<VoronoiMethods<Type, PType>>  vmethods;

static const int DEGREE = 3;

PType root_list[3];
PType rho1, rho2, max_minima_dist; // 2nd and 3rd elements of the integral basis

Type discriminant;
Type E =Type(0); Type G = Type(0);

long index; // The number f such that disc(O) = f^2 * Delta, Delta = field discriminant
int unit_rank, r1,r2;

// See the issue link below regarding the four members below
// https://gitlab.cpsc.ucalgary.ca/jacobs/ANTL/-/issues/11
Type class_number = Type(0);
std::vector<Type> cg_structure;

PType regulator = PType(0);
std::vector<CubicElement<Type, PType>> fundamentalUnits;

//std::vector<std::vector<PType>(DEGREE, PType(0)) > embedding_matrix;

Type mul_table[3][3];

static PType order_temp;
CubicElement<Type, PType> temp_element = CubicElement<Type, PType>(this, Type(1), Type(0), Type(0), Type(1));


void set_roots();

// do I even need this?
void set_mul_table();

void set_integral_basis();

void set_class_number();

void set_class_group();

void set_unit_rank();

void set_max_minima_dist();

virtual void set_regulator() = 0;

//void compute_fundamental_unit();

private:

CubicOrder();   // private default constructor means you can't call this




}; // closes class scope


template<typename T, typename PT> PT CubicOrder<T,PT>::order_temp;
#include "../../../src/Cubic/CubicOrder.cpp"

#endif // guard
