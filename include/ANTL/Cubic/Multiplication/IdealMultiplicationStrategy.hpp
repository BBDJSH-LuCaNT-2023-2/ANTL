#ifndef ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP
#define ANTL_IDEAL_MULTIPLICATION_STRATEGY_HPP

// strategy class for ideal multiplications


template<typename Type, typename PType>
class IdealMultiplicationStrategy {

public:


    // this method should do some checking to ensure that the ideals are in the same order

    //this will be a virtual function which is instantiated in subclasses
    // different subclasses shall implement distinct methods

//virtual void multiply(CubicIdeal<Type> &C, const CubicIdeal<Type> &A, const CubicIdeal<Type> &B) = 0;



protected:

// variables to hold the product of two ideals as an integral 3 by 3 matrix along with a denominator.
Type product_basis[3][3];
Type denom;

private:



}; //end class def

#endif // include guard
