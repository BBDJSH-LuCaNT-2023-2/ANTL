#ifndef ANTL_CUBIC_ORDER_REAL_CPP
#define ANTL_CUBIC_ORDER_REAL_CPP

#include "../../include/ANTL/Cubic/CubicOrderReal.hpp"



template <typename Type, typename PType>
void CubicOrderReal<Type, PType> :: compute_conjugate_elements(){

  mul(this->conjugate_bases[0][0],to<PType>((this->defining_IBCF)[3]), this->root_list[1]);           // rho1' = a*delta'

  this->conjugate_bases[0][1] = this->conjugate_bases[0][0];                                          // a*delta'
  add(this->conjugate_bases[0][1],this->conjugate_bases[0][1],to<PType>((this->defining_IBCF)[2]) );  // a*delta' + b
  mul(this->conjugate_bases[0][1], this->conjugate_bases[0][1], this->root_list[1]);                  // rho2' = a*delta'^2 + b*delta'


  mul(this->conjugate_bases[1][0],to<PType>((this->defining_IBCF)[3]), this->root_list[2]);           // a*delta''
  this->conjugate_bases[1][1] = this->conjugate_bases[1][0];                                          // a*delta''
  add(this->conjugate_bases[1][1],this->conjugate_bases[1][1], to<PType>( (this->defining_IBCF)[2]) ); // a*delta'' + b
  mul(this->conjugate_bases[1][1], this->conjugate_bases[1][1], this->root_list[2]);                   // a*delta''^2 + b*delta''

  #ifdef DEBUG2
  std::cout << "Computing conjugate elements" << this->get_conjugate_bases(0,0) << " " << this->get_conjugate_bases(0,1)<< std::endl;
  std::cout << "Computing conjugate elements" << this->get_conjugate_bases(1,0) << " " << this->get_conjugate_bases(1,1)<< std::endl;
  #endif
}

template <typename Type, typename PType>
void CubicOrderReal<Type, PType> :: get_real_value(PType & newVal, const Type &U, const Type &X, const Type &Y, const Type &D, int conj){
    //std::cout << "real order, real value" << std::endl;
    NTL::clear(newVal);
    if (this->is_real() && (conj != 0)){
      //std::cout << this->get_conjugate_bases(conj,0) << std::endl;
       --conj;
       //std::cout << "trying to get conjugate base elements" << std::endl;
       //std::cout << this->get_conjugate_bases(conj,1) << " " << this->get_conjugate_bases(conj,0) << std::endl;

       mul(newVal, to<PType>(X), this->get_conjugate_bases(conj,0));  // x*rho1'

       add(newVal, newVal, to<PType>(U));

       mul(this->order_temp, to<PType>(Y), this->get_conjugate_bases(conj,1)); // y* rho2'

       add(newVal, newVal, this->order_temp);
    }
    else{

      mul(newVal, to<PType>(X), this->get_rho1());
      add(newVal, newVal, to<PType>(U));
      mul(this->order_temp, to<PType>(Y), this->get_rho2());
      add(newVal, newVal, this->order_temp);
    }
    div(newVal, newVal, to<PType>(D));

}


template <typename Type, typename PType>
void CubicOrderReal<Type, PType> ::set_regulator(){
  if (this->regulator == 0){
    if(this->fundamentalUnits.size() == 0){
        this->unit_strat->compute(this->fundamentalUnits, this, this->is_real());
    }
    PType reg_temp, reg_temp2;

    //
    //std::cout << "regulator calc " << this->fundamentalUnits[0].get_u() << this->fundamentalUnits[0].get_x() <<this->fundamentalUnits[0].get_y() <<this->fundamentalUnits[0].get_denom()<< std::endl;
    this->get_real_value(reg_temp2, this->fundamentalUnits[0].get_u(),this->fundamentalUnits[0].get_x(),this->fundamentalUnits[0].get_y(),this->fundamentalUnits[0].get_denom(), 0);
    abs(reg_temp2, reg_temp2);
    log(reg_temp2, reg_temp2);                            // log( | sigma_0(e1) |)
    this->get_real_value(reg_temp, this->fundamentalUnits[1].get_u(),this->fundamentalUnits[1].get_x(),this->fundamentalUnits[1].get_y(), this->fundamentalUnits[1].get_denom(), 1);
    abs(reg_temp, reg_temp);
    log(reg_temp, reg_temp);                                             // log( | sigma_1(e2) |)
    mul(this->regulator, reg_temp2, reg_temp);

    this->get_real_value(reg_temp2, this->fundamentalUnits[0].get_u(),this->fundamentalUnits[0].get_x(),this->fundamentalUnits[0].get_y(), this->fundamentalUnits[0].get_denom(), 1);
    abs(reg_temp2, reg_temp2);
    log(reg_temp2, reg_temp2);                            // log( | sigma_1(e1) |)
    this->get_real_value(reg_temp, this->fundamentalUnits[1].get_u(),this->fundamentalUnits[1].get_x(),this->fundamentalUnits[1].get_y(),this->fundamentalUnits[1].get_denom(), 0);
    abs(reg_temp, reg_temp);
    log(reg_temp, reg_temp);                                            // log( | sigma_0(e2) |)
    mul(reg_temp2, reg_temp2, reg_temp);

    sub(this->regulator, this->regulator, reg_temp2);            // determinant of log matrix
    abs(this->regulator, this->regulator);

  }
}






#endif
