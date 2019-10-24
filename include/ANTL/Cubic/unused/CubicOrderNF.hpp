#ifndef ANTL_CUBIC_ORDER_NF_H
#define ANTL_CUBIC_ORDER_NF_H

#include "CubicOrder.hpp"
//#include "Multiplication/IdealMultiplicationStrategy.hpp"


using boost::math::tools::polynomial;
using namespace ANTL;


// forward declaration
template<typename Type, typename PType>
class CubicOrder;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class CubicIdeal;
template<typename Type, typename PType>
class IdealMultiplicationStrategy;


template <typename Type, typename PType>
class CubicOrderNF : public CubicOrder<Type, PType>{

public:



// Class methods

/**
* @brief Constructor which requires a defining polynomial as well as a multiplication strategy
* @param A cubic polynomial, a pointer to a IdealMultiplicationStrategy object
*/
//CubicOrderNF<Type,PType> (polynomial<Type> const &poly, IdealMultiplicationStrategy<Type, PType> * mulstrat);
CubicOrderNF<Type,PType> (polynomial<Type> const &poly);


bool is_equal(const CubicOrder<Type, PType> &CO2) const;

void roots_swap_position(int p1, int p2);

int splitting_type(int p);

/**
* @brief Method for setting the IdealMultiplicationStrategy
* @param A pointer to a concrete strategy (subclass of IdealMultiplicationStrategy)
*/
inline void set_mul_strategy(IdealMultiplicationStrategy<Type,PType> * mstrat){
  this->mul_method = mstrat;
};


void mul(CubicIdeal<Type,PType> & A, const CubicIdeal<Type,PType> & B, const CubicIdeal<Type,PType> & C);


/********************** friend functions **********************/
template <typename T, typename PT>
friend bool is_equal(const CubicOrderNF<T, PT> &CO1, const CubicOrderNF<T, PT> &CO2);


protected:

// Class variables    REQUIRE(ord2.is_equal(Rufio) == true);


// Each column of mulTable represents the product of two canonical basis elements in terms of the
// canonical integral basis itself.
// The order is rho1^2, rho1*rho2, rho2^2

Type mul_table[3][3];

std::vector<PType> fundamentalUnits;

void set_roots();


// do I even need this?
void set_mul_table();

void set_integral_basis();

void set_class_number();

void set_class_group();

void set_regulator();

private:




}; //close class definition
#include "../../../src/Cubic/CubicOrderNF.cpp"
#endif // guard
