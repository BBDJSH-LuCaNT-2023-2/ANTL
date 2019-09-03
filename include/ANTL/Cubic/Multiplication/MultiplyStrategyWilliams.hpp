#ifndef ANTL_MULTIPLY_WILLIAMS_HPP
#define ANTL_MULTIPLY_WILLIAMS_HPP


template<typename Type, typename PType>
class MultiplyStrategyWilliams : public IdealMultiplicationStrategy<Type, PType> {

public:

void multiply(CubicIdeal<Type> &C, const CubicIdeal<Type> &A, const CubicIdeal<Type> &B){
  std::cout <<"Williams multiply" << std::endl;
}



protected:





private:

Type AMatrix1;
Type AMatrix2;

productBasis = Type[3][3];  //place to store product coefficients
xVector = Type[9]; // holds the Bezout coefficients of C_3
genSet: Type[3][9];   //This is unavoidable since C_3 is the gcd(c_1, ... c_9)
                      // and A_3 and B_3 are defined in terms of the a_i, b_i and the above gcd bezout coefficients.

// at this point I think it is possible to determine the product basis without
// introducing any further variables.




}




#endif
