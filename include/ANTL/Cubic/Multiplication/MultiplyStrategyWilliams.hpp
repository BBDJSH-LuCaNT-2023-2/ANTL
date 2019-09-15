#ifndef ANTL_MULTIPLY_WILLIAMS_HPP
#define ANTL_MULTIPLY_WILLIAMS_HPP

#include "IdealMultiplicationStrategy.hpp"
#include "../CubicIdealNF.hpp"
//template<typename Type, typename PType>
//class CubicIdealNF;

template<typename Type, typename PType>
class MultiplyStrategyWilliams : public IdealMultiplicationStrategy<Type, PType> {

public:

  MultiplyStrategyWilliams();

  virtual ~MultiplyStrategyWilliams(){
    //delete mul_holder;
    delete this->my_order;
  };

void multiply(CubicIdealNF<Type,PType> &A, const CubicIdealNF<Type,PType> &B, const CubicIdealNF<Type,PType> &C);


protected:





private:

CubicElementNF<Type, PType> mul_holder = CubicElementNF<Type, PType>();

Type bezout1;
Type product_basis[3][3];  //place to store product coefficients
Type x_vector[9]; // holds the Bezout coefficients of C_3
Type gen_set[3][9];   //This is unavoidable since C_3 is the gcd(c_1, ... c_9)
                      // and A_3 and B_3 are defined in terms of the a_i, b_i and the above gcd bezout coefficients.

Type c_set[9];
// at this point I think it is possible to determine the product basis without
// introducing any further variables.




}; //end class definition
#include "../../../../src/Cubic/Multiplication/MultiplyStrategyWilliams.cpp"



#endif
