#ifndef GUARD_ztest_cpp
#define GUARD_ztest_cpp

#include <iostream>

// qvm matrix headers
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>

#include "../../../include/ANTL/Cubic/generalFunctions.hpp"
#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/CubicOrderNF.hpp"
#include "../../../include/ANTL/Cubic/CubicElementNF.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
using namespace NTL;
using namespace ANTL;
using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using std::cout;
using std::endl;
// For floating point arithmetic error tolerance
const double DOUBLE_TOLERANCE = 0.0000001;

using namespace NTL;

NTL_CLIENT
//using namespace NTL;
// 020-TestCase.cpp
// Let Catch provide main():
//#include "../catch.hpp"

int main(){

  ZZ joe, bob;
  joe = 1;
  bob = 2;

  std::cout <<joe+bob << std::endl;

}

#endif
