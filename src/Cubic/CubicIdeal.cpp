#ifndef ANTL_CUBIC_IDEAL_CPP
#define ANTL_CUBIC_IDEAL_CPP

#include "../../include/ANTL/Cubic/CubicIdeal.hpp"


template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: normalize(){
  //set the new denominator
  mul(this->ci_temp, gen1.get_denominator(), gen2.get_denominator());
  mul(this->ci_temp, this->ci_temp, gen3.get_denominator()); //ci_temp = d1*d2*d3
  mul(this->denom, this->denom, this->ci_temp);


  // Note we should be careful here since mul automatically normalizes after
  // completing the operation. By multiplying by d1*d2*d3, the denominators
  // of each generator should be 1. so there should be no problems normalizing
  // ideal
  mul(gen1, this->ci_temp);
  mul(gen2, this->ci_temp);
  mul(gen3, this->ci_temp);

  this->ci_temp = GCD(GCD(gen1.get_u(), gen2.get_u()),gen3.get_u() );
  if (ci_temp != Type(1)){
    this->ci_temp = GCD(GCD(GCD(GCD(GCD(GCD(GCD(
      this->ci_temp,
      gen1.get_x() ),
      gen2.get_x() ),
      gen3.get_x() ),
      gen1.get_y() ),
      gen2.get_y() ),
      gen3.get_y() ),
      this->denom  );

      div(gen1, gen1, ci_temp);
      div(gen2, gen2, ci_temp);
      div(gen3, gen3, ci_temp);

      div(this->denom, this->denom, this->ci_temp);
  }

}

#endif
