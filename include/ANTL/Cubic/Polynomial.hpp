#ifndef ANTL_POLYNOMIAL_HPP
#define ANTL_POLYNOMIAL_HPP


#include <NTL/ZZX.h>
//using NTL::ZZX; using NTL::SetCoeff;
using namespace NTL;

// It is intended that Type is the datatype for integers
// and that PType holds real numbers.
template <typename Type, typename PType>
class Polynomial {

public:
  // constructors
  Polynomial(Type A, Type B, Type C, Type D){
      coefficients.assign(A,B,C,D);
  };






private:
  vector<Type> coefficients;


protected:


};


#endif // guard
