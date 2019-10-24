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

  }


void roots_swap_position(int p1, int p2){
  if (this->discriminant > 0) {
    // swap roots and recalculate the integral basis.
    std::swap(this->root_list[p1], this->root_list[p2]);
    this->set_integral_basis();
    compute_conjugate_elements();
  }

};


void compute_conjugate_elements();

inline PType get_conjugate_bases(int i, int j) const {
  return conjugate_bases[i][j];
}

void get_real_value(PType & newVal, Type &U, Type &X, Type &Y, Type &D, int conj = 0);


/**
* @brief See pg 317 of CFG, this algorithm computes a pair of fund. units of a cubic order
*/
void compute_fundamental_unit();

CubicElement<Type, PType> * get_fund_unit(int i){
  if (i > 1){
    std::cout << "Please supply argument either 0 or 1" << std::endl;
  }
  if (this->fundamentalUnits.size() == 0){
    this->compute_fundamental_unit();
    return &this->fundamentalUnits[i];
  }
  return &this->fundamentalUnits[i];
};


protected:
// This is meant to store the Matrix
// [ rho1', rho2' ]
// [ rho1'' rho2'']
PType conjugate_bases[2][2];
std::vector<CubicIdeal<Type, PType>> x_cycle;         // container to hold x-cycle
std::vector<CubicElement<Type, PType>> adj_minima_vec;    // container for holding adjacent min


private:



};



#include "../../../src/Cubic/CubicOrderReal.cpp"

#endif
