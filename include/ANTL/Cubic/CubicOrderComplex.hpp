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

  PType interim = to<PType>(0);

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

// this should be analogous to the NEAREST or TARGET algorithm of Thiel, or Reduce algorithm of Schoof
// move this to th cpp file later
void close_minimum(CubicIdeal<Type, PType> & reduced_ideal, CubicElement<Type, PType> & minimum, CubicIdeal<Type, PType> & ideal1, const std::vector<PType> & vec1){
  // idea: vec1 should be a real type vector of length r+1
  // want to find a minimum in ideal1 which is close to vec1, so use voronoi


  PType minima_log, step_log, abs_log, debug_var;
  CubicElement<Type, PType> step_element = CubicElement<Type, PType>(this);
  ideal1.reduce(reduced_ideal, minimum);
  cout << "\nClose Starts" << endl;
  minimum.get_real_value(minima_log);

  log(minima_log, minima_log);
  cout << "minima_log" << minima_log << endl;
  abs(abs_log, minima_log);

  if(abs_log > vec1[0]){
    cout << "Minimum is too big!" << endl;
    return;
  }else{
    reduced_ideal.make_voronoi_basis();
    step_element.assign(reduced_ideal.get_coeff(0,1),reduced_ideal.get_coeff(1,1),reduced_ideal.get_coeff(2,1), reduced_ideal.get_coeff(0,0));
    step_element.get_real_value(step_log);
    log(step_log, step_log);
    add(step_log, step_log, minima_log);
    abs(abs_log, step_log);

    while(abs_log < vec1[0]){
        reduced_ideal.divide_adjacent(reduced_ideal, step_element);
        reduced_ideal.make_canonical();
        minima_log = step_log;

        mul(minimum, minimum, step_element);

        minimum.get_real_value(debug_var);
        log(debug_var, debug_var);
        cout << "log of the minimum: " << debug_var << endl;
        cout << "minima_log " << step_log << endl;

        reduced_ideal.make_voronoi_basis();
        step_element.assign(reduced_ideal.get_coeff(0,1),reduced_ideal.get_coeff(1,1),reduced_ideal.get_coeff(2,1), reduced_ideal.get_coeff(0,0));
        step_element.get_real_value(step_log);
        log(step_log, step_log);
        add(step_log, step_log, minima_log);
        abs(abs_log, step_log);

        cout << "How close are we? " << step_log << "  " << vec1[0] << endl;


    }
    cout << "Final: " << minima_log << "  " << vec1[0] << endl;
  }

};
protected:




private:



};


//#include "../../../src/Cubic/CubicOrderComplex.cpp"
#endif // guard
