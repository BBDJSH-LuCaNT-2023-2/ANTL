#ifndef ANTL_MULTIPLY_WILLIAMS_CPP
#define ANTL_MULTIPLY_WILLIAMS_CPP

#include "../../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"


template<typename Type, typename PType>
MultiplyStrategyWilliams<Type, PType> :: MultiplyStrategyWilliams(){

}


template<typename Type, typename PType>
void MultiplyStrategyWilliams<Type, PType> :: multiply(CubicIdealNF<Type,PType>  &A, const CubicIdealNF<Type,PType>  &B, const CubicIdealNF<Type,PType> &C){
  mul_holder.set_order(A.get_order());

  mul(this->mul_holder, *A.get_gen1(),*B.get_gen1());
  this->gen_set[0][0] = this->mul_holder.get_u();
  this->gen_set[1][0] = this->mul_holder.get_x();
  this->gen_set[2][0] = this->mul_holder.get_y();

  mul(this->mul_holder, *A.get_gen1(),*B.get_gen2());
  this->gen_set[0][1] = this->mul_holder.get_u();
  this->gen_set[1][1] = this->mul_holder.get_x();
  this->gen_set[2][1] = this->mul_holder.get_y();

  mul(this->mul_holder, *A.get_gen1(),*B.get_gen3());
  this->gen_set[0][2] = this->mul_holder.get_u();
  this->gen_set[1][2] = this->mul_holder.get_x();
  this->gen_set[2][2] = this->mul_holder.get_y();

  mul(this->mul_holder, *A.get_gen2(),*B.get_gen1());
  this->gen_set[0][3] = this->mul_holder.get_u();
  this->gen_set[1][3] = this->mul_holder.get_x();
  this->gen_set[2][3] = this->mul_holder.get_y();
  mul(this->mul_holder, *A.get_gen2(),*B.get_gen2());
  this->gen_set[0][4] = this->mul_holder.get_u();
  this->gen_set[1][4] = this->mul_holder.get_x();
  this->gen_set[2][4] = this->mul_holder.get_y();
  mul(this->mul_holder, *A.get_gen2(),*B.get_gen3());
  this->gen_set[0][5] = this->mul_holder.get_u();
  this->gen_set[1][5] = this->mul_holder.get_x();
  this->gen_set[2][5] = this->mul_holder.get_y();


  mul(this->mul_holder, *A.get_gen3(),*B.get_gen1());
  this->gen_set[0][6] = this->mul_holder.get_u();
  this->gen_set[1][6] = this->mul_holder.get_x();
  this->gen_set[2][6] = this->mul_holder.get_y();
  mul(this->mul_holder, *A.get_gen3(),*B.get_gen2());
  this->gen_set[0][7] = this->mul_holder.get_u();
  this->gen_set[1][7] = this->mul_holder.get_x();
  this->gen_set[2][7] = this->mul_holder.get_y();
  mul(this->mul_holder, *A.get_gen3(),*B.get_gen3());
  this->gen_set[0][8] = this->mul_holder.get_u();
  this->gen_set[1][8] = this->mul_holder.get_x();
  this->gen_set[2][8] = this->mul_holder.get_y();

  mul(this->denom, A.get_denom(), B.get_denom());

  // ****************** find C_3, the gcd of the bottom row ***************** //

  //get the gcd of the first two entries, and place them in the lower right corner
  XGCD(this->product_basis[2][2], x_vector[0], x_vector[1], this->gen_set[2][0], this->gen_set[2][1]);

  // each outer loop iteration finds the gcd of the previous gcd and the next c_i
  // In other words, g_1 = gcd(c_0,c_1), g_i = gcd(g_(i-1), c_i) for 9 > i >= 2
  // The bezout coefficients x,y such that x*g_{i-1} + y*c_i = g_i must be tracked
  // y is placed in the array, while x is used to update the previous coefficients
  for (int i= 2; i < 9; ++i){
    XGCD(this->product_basis[2][2], bezout1 , x_vector[i], this->product_basis[2][2], this->gen_set[2][i]);

    // update previous coefficients
    for (int j = 0; j<i; ++j){
      x_vector[i]*= bezout1;
    }
  }
  long check=0;
  for (int i = 0; i < 9; i++){
    check += x_vector[i]*this->gen_set[2][i];
  }
  std::cout << check << " vs " << this->product_basis[2][2] << std::endl;


} // close function mul








#endif
