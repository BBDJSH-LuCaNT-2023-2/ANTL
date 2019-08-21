#ifndef ANTL_CUBIC_ORDER_NF_H
#define ANTL_CUBIC_ORDER_NF_H


#include <NTL/ZZX.h>
#include <boost/math/tools/polynomial.hpp>

#include "CubicOrder.hpp"
#include "GeneralTemplateFunctions.hpp"
#include "Multiplication/IdealMultiplicationStrategy.hpp"
using boost::math::tools::polynomial;


template <typename Type, typename PType>
class CubicOrderNF : public CubicOrder<Type, PType>{

public:



// Class methods

// constructor
CubicOrderNF<Type,PType> (polynomial<Type> const &poly);



bool is_equal(const CubicOrder<Type, PType> &CO2) const;


int splitting_type(int p);

/********************** friend functions **********************/
template <typename T, typename PT>
friend bool is_equal(const CubicOrder<T, PT> &CO1, const CubicOrder<T, PT> &CO2);


protected:

// Class variables


// Each column of mulTable represents the product of two canonical basis elements in terms of the
// canonical integral basis itself.
// The order is rho1^2, rho1*rho2, rho2^2

Type mul_table[3][3];

std::vector<PType> fundamentalUnits;

IdealMultiplicationStrategy<Type, PType> * mul_method;

void set_mul_table();

void set_integral_basis();

void set_class_number();

void set_class_group();

void set_regulator();
private:




}; //close class definition
#include "../../../src/Cubic/CubicOrderNF.cpp"
#endif // guard
