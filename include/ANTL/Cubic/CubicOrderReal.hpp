#ifndef ANTL_CUBIC_ORDER_REAL_H
#define ANTL_CUBIC_ORDER_REAL_H


#include "CubicOrder.hpp"

template<typename Type, typename PType>
class CubicOrder;

template<typename Type, typename PType>
class CubicIdeal;

template<typename Type, typename PType>
class CubicElement;

template <typename T, typename PT>
void mul (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C);
template <typename T, typename PT>
bool is_equal(CubicIdeal <T,PT> & A, CubicIdeal <T,PT> & B);

template <typename Type, typename PType>
class CubicOrderReal : public CubicOrder<Type, PType> {

public:


CubicOrderReal( polynomial<Type> const &poly)
  : CubicOrder<Type,PType>::CubicOrder(poly) {

    this->set_integral_basis();
    compute_conjugate_elements();
  }

// Swaps the roots as well as the integral basis values with with specified conjugate
void roots_swap_position(int p1, int p2){
  if (this->discriminant > 0) {
    // swap roots and recalculate the integral basis.
    std::swap(this->root_list[p1], this->root_list[p2]);
    std::swap(this->rho1, this->conjugate_bases[p2-1][0]);
    std::swap(this->rho2, this->conjugate_bases[p2-1][1]);
    this->set_integral_basis();
    compute_conjugate_elements();
  }

};


void compute_conjugate_elements();

inline PType get_conjugate_bases(int i, int j) const {
  return conjugate_bases[i][j];
}

void get_real_value(PType & newVal, const Type &U, const Type &X, const Type &Y, const Type &D, int conj = 0);


/**
* @brief See pg 317 of CFG, this algorithm computes a pair of fund. units of a cubic order
*/
CubicElement<Type, PType> * get_fundamental_unit(int i){
  if (this->fundamentalUnits.size() == 0){
    this->unit_strat->compute(this->fundamentalUnits, this, this->is_real());
    return &this->fundamentalUnits[i];
  }else
  return &this->fundamentalUnits[i];
};


void set_regulator();

protected:
// This is meant to store the Matrix
// [ rho1', rho2' ]
// [ rho1'' rho2'']
PType conjugate_bases[2][2];



private:



};



#include "../../../src/Cubic/CubicOrderReal.cpp"

#endif
