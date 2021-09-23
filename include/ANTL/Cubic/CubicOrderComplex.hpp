#ifndef ANTL_CUBIC_ORDER_COMPLEX_H
#define ANTL_CUBIC_ORDER_COMPLEX_H

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
class CubicOrderComplex : public CubicOrder<Type, PType> {

public:


CubicOrderComplex( polynomial<Type> const &poly)
  : CubicOrder<Type,PType>::CubicOrder(poly) {

  }


void get_real_value(PType & newVal, const Type &U, const Type &X, const Type &Y, const Type &D, int conj = 0){

  PType interim = to<PType>(0);         // formerly used static variable this->order_temp

  NTL::mul(interim, to<PType>(Y), this->get_rho2());
  NTL::mul(newVal, to<PType>(X), this->get_rho1());
  NTL::add(newVal, newVal, to<PType>(U));
  NTL::add(newVal, interim, newVal);
  NTL::div(newVal, newVal, to<PType>(D));
}


/**
* @brief See pg 317 of CFG, this algorithm returns the specified fund. unit of a cubic order
*/
CubicElement<Type, PType> * get_fundamental_unit(int i){

if (this->fundamentalUnits.size() == 0){
  this->unit_strat->compute(this->fundamentalUnits, this, this->is_real());
  return &this->fundamentalUnits[0];
}else{
  return &this->fundamentalUnits[0];
}

};


void set_regulator(){

  if (this->regulator == 0){
    if(this->fundamentalUnits.size() == 0){
        this->unit_strat->compute(this->fundamentalUnits, this, this->is_real());
    }
    get_real_value(this->order_temp, this->fundamentalUnits[0].get_u(),this->fundamentalUnits[0].get_x(),this->fundamentalUnits[0].get_y(), this->fundamentalUnits[0].get_denom() );
    abs(this->order_temp, this->order_temp);
    log(this->order_temp, this->order_temp);
    this->regulator = this->order_temp;
  }
}

protected:




private:



};


//#include "../../../src/Cubic/CubicOrderComplex.cpp"
#endif // guard
