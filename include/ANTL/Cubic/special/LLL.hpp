#ifndef CUBIC_LLL_HPP
#define CUBIC_LLL_HPP


#include <boost/math/bindings/rr.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/math/tools/polynomial.hpp>

#include "../../Arithmetic/QQ.hpp"
#include "../../common.hpp"

#include <functional>
#include <cstdint>
#include <iostream>
#include <NTL/RR.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
using namespace NTL;
using namespace ANTL;


template<PType> void LLL_reduction(PType logbasis[2][2]);

template<PType> void LLL_reduction(PType lambda00, PType lambda10, PType lambda01, PType lambda11);
//const quad_float eps=1e-23;
//const quad_float PI =    3.14159265358979323846264;
//const	quad_float TwoPi = 6.28318530717958647692528;






#include "../../../src/Cubic/special/LLL.cpp"


#endif
