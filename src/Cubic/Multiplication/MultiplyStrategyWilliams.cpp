#ifndef ANTL_MULTIPLY_WILLIAMS_CPP
#define ANTL_MULTIPLY_WILLIAMS_CPP

#include "../../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"


template<typename Type, typename PType>
MultiplyStrategyWilliams<Type, PType> :: MultiplyStrategyWilliams(){

}


template<typename Type, typename PType>
void MultiplyStrategyWilliams<Type, PType> :: multiply(CubicIdealNF<Type,PType>  &A, const CubicIdealNF<Type,PType>  &B, const CubicIdealNF<Type,PType> &C){
  mul_holder.set_order(A.get_order());


  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ****************** Compute the set of 9 generators ********************* //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  mul(this->mul_holder, *B.get_gen1(),*C.get_gen1());
  this->gen_set[0][0] = this->mul_holder.get_u();
  this->gen_set[1][0] = this->mul_holder.get_x();
  this->gen_set[2][0] = this->mul_holder.get_y();

  mul(this->mul_holder, *B.get_gen1(),*C.get_gen2());
  this->gen_set[0][1] = this->mul_holder.get_u();
  this->gen_set[1][1] = this->mul_holder.get_x();
  this->gen_set[2][1] = this->mul_holder.get_y();

  mul(this->mul_holder, *B.get_gen1(),*C.get_gen3());
  this->gen_set[0][2] = this->mul_holder.get_u();
  this->gen_set[1][2] = this->mul_holder.get_x();
  this->gen_set[2][2] = this->mul_holder.get_y();

  mul(this->mul_holder, *B.get_gen2(),*C.get_gen1());
  this->gen_set[0][3] = this->mul_holder.get_u();
  this->gen_set[1][3] = this->mul_holder.get_x();
  this->gen_set[2][3] = this->mul_holder.get_y();
  mul(this->mul_holder, *B.get_gen2(),*C.get_gen2());
  this->gen_set[0][4] = this->mul_holder.get_u();
  this->gen_set[1][4] = this->mul_holder.get_x();
  this->gen_set[2][4] = this->mul_holder.get_y();
  mul(this->mul_holder, *B.get_gen2(),*C.get_gen3());
  this->gen_set[0][5] = this->mul_holder.get_u();
  this->gen_set[1][5] = this->mul_holder.get_x();
  this->gen_set[2][5] = this->mul_holder.get_y();


  mul(this->mul_holder, *B.get_gen3(),*C.get_gen1());
  this->gen_set[0][6] = this->mul_holder.get_u();
  this->gen_set[1][6] = this->mul_holder.get_x();
  this->gen_set[2][6] = this->mul_holder.get_y();
  mul(this->mul_holder, *B.get_gen3(),*C.get_gen2());
  this->gen_set[0][7] = this->mul_holder.get_u();
  this->gen_set[1][7] = this->mul_holder.get_x();
  this->gen_set[2][7] = this->mul_holder.get_y();
  mul(this->mul_holder, *B.get_gen3(),*C.get_gen3());
  this->gen_set[0][8] = this->mul_holder.get_u();
  this->gen_set[1][8] = this->mul_holder.get_x();
  this->gen_set[2][8] = this->mul_holder.get_y();

  mul(this->denom, B.get_denom(), C.get_denom());


  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ****************** Find C_3, the gcd of the bottom row ***************** //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  //get the gcd of the first two entries, and place them in the lower right corner
  XGCD(this->product_basis[2][2], x_vector[0], x_vector[1], this->gen_set[2][0], this->gen_set[2][1]);

  // Outer loop iteration finds gcd of  previous gcd and the next c_i
  // In other words, g_1 = gcd(c_0,c_1), g_i = gcd(g_(i-1), c_i) for 9 > i >= 2

  // The bezout coefficients x,y such that x*g_{i-1} + y*c_i = g_i must be tracked
  // y is placed in the array, while x is used to update the previous coefficients
  for (int i= 2; i < 9; ++i){
    XGCD(this->product_basis[2][2], bezout1 , x_vector[i], this->product_basis[2][2], this->gen_set[2][i]);

    // update previous coefficients
    for (int j = 0; j<i; ++j){
      mul(x_vector[i], x_vector[i], bezout1);
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ******************         Compute A_3 and B_3     ********************* //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////

  ::clear(product_basis[0][2]);
  ::clear(product_basis[1][2]);

  for (int i = 0; i < 9; ++i){
    mul(bezout1, x_vector[i], gen_set[1][i]);
    add(this->product_basis[1][2], this->product_basis[1][2], bezout1); //computes B_3

    mul(bezout1, x_vector[i], gen_set[0][i]);
    add(this->product_basis[0][2], this->product_basis[0][2], bezout1); // computes A_3

  }

  // updating the entries
  // a -> a' = a_i - (A_3 * c_i / C_3)
  // b -> b' = b_i - (B_3 * c_i / C_3)
  for (int i = 0; i < 9; ++i){
    div(bezout1,this->gen_set[2][i], this->product_basis[2][2]); // store c_i/C3 in bezout1

    // note that product_basis[0][0] is being used as an intermediate variable holder here
    mul(this->product_basis[0][0], bezout1, this->product_basis[1][2]); // B_3 * c_i / C_3
    sub(gen_set[1][i], gen_set[1][i], this->product_basis[0][0]); // b_i - B_3 * c_i / C_3

    mul(this->product_basis[0][0], bezout1, this->product_basis[0][2]); // A_3 * c_i / C_3
    sub(gen_set[0][i], gen_set[0][i], this->product_basis[0][0]); // a_i - A_3 * c_i / C_3


  }

  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ****************** Compute B_2 and A_2 (second basis vector) *********** //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  ::clear(this->product_basis[0][1]);
  ::clear(this->product_basis[1][1]);
  ::clear(this->product_basis[2][1]);

  XGCD(this->product_basis[1][1], x_vector[0], x_vector[1], this->gen_set[1][0], this->gen_set[1][1]);

  for (int i= 2; i < 9; ++i){
    XGCD(this->product_basis[1][1], bezout1 , x_vector[i], this->product_basis[1][1], this->gen_set[1][i]);

    // update previous coefficients
    for (int j = 0; j<i; ++j){
      mul(x_vector[i], x_vector[i], bezout1);
    }
  }

  for (int i = 0; i < 9; ++i){
    mul(bezout1, x_vector[i], gen_set[0][i]);
    add(this->product_basis[0][1], this->product_basis[0][1], bezout1); // computes A_2

  }

  // updating the entries
  // a' -> a'' = a'_i - (A_2 * b'_i / B_2)
  for (int i = 0; i < 9; ++i){
    div(bezout1,this->gen_set[1][i], this->product_basis[1][1]); // store b'_i / B2 in bezout1

    // note that product_basis[0][0] is being used as an intermediate variable holder here

    mul(this->product_basis[0][0], bezout1, this->product_basis[0][1]); // A_2 * b'_i / B_2
    sub(gen_set[0][i], gen_set[0][i], this->product_basis[0][0]); // a'_i - (A_2 * b'_i / B_2)


  }
  //////////////////////////////////////////////////////////////////////////////
  //                                                                          //
  // ****************** Compute A_3 (the first basis vector) **************** //
  //                                                                          //
  //////////////////////////////////////////////////////////////////////////////
  ::clear(this->product_basis[0][0]);
  ::clear(this->product_basis[1][0]);
  ::clear(this->product_basis[2][0]);


  XGCD(this->product_basis[0][0],x_vector[0],x_vector[1], this->gen_set[0][0], this->gen_set[0][1]);

  for (int i= 2; i < 9; ++i){
    XGCD(this->product_basis[0][0], bezout1,x_vector[i], this->product_basis[0][0], this->gen_set[0][i]);

  }

  A.assign(this->product_basis[0][0],this->product_basis[1][0],this->product_basis[2][0], this->denom,
  this->product_basis[0][1],this->product_basis[1][1],this->product_basis[2][1], this->denom,
  this->product_basis[0][2],this->product_basis[1][2],this->product_basis[2][2], this->denom);
} // close function mul








#endif
