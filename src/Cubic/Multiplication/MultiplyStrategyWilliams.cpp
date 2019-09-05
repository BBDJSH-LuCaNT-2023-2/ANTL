#ifndef ANTL_MULTIPLY_WILLIAMS_CPP
#define ANTL_MULTIPLY_WILLIAMS_CPP

#include "../../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"


template<typename Type, typename PType>
void MultiplyStrategyWilliams<Type, PType> :: mul(CubicIdeal<Type,PType>  &A, const CubicIdeal<Type,PType>  &B, const CubicIdeal<Type,PType> &C){

  mul(this->mul_holder, A.gen1,B.gen1);
  this->gen_set[0][0] = this->mul_holder.get_u();
  this->gen_set[1][0] = this->mul_holder.get_x();
  this->gen_set[2][0] = this->mul_holder.get_y();

  mul(this->mul_holder, A.gen1,B.gen2);
  this->gen_set[0][1] = this->mul_holder.get_u();
  this->gen_set[1][1] = this->mul_holder.get_x();
  this->gen_set[2][1] = this->mul_holder.get_y();

  mul(this->mul_holder, A.gen1,B.gen3);
  this->gen_set[0][2] = this->mul_holder.get_u();
  this->gen_set[1][2] = this->mul_holder.get_x();
  this->gen_set[2][2] = this->mul_holder.get_y();

  mul(this->mul_holder, A.gen2,B.gen2);
  this->gen_set[0][3] = this->mul_holder.get_u();
  this->gen_set[1][3] = this->mul_holder.get_x();
  this->gen_set[2][3] = this->mul_holder.get_y();
  mul(this->mul_holder, A.gen1,B.gen3);
  this->gen_set[0][4] = this->mul_holder.get_u();
  this->gen_set[1][4] = this->mul_holder.get_x();
  this->gen_set[2][4] = this->mul_holder.get_y();
  mul(this->mul_holder, A.gen3,B.gen3);
  this->gen_set[0][5] = this->mul_holder.get_u();
  this->gen_set[1][5] = this->mul_holder.get_x();
  this->gen_set[2][5] = this->mul_holder.get_y();

  mul(this->denom, A.get_denom(), B.get_denom());
}








#endif
