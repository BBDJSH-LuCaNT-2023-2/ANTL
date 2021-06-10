#ifndef ANTL_CUBIC_HEADER_H
#define ANTL_CUBIC_HEADER_H

#include <iostream>
#include <functional>
#include <iterator>
#include <unordered_map>
// qvm matrix headers
//#include <boost/qvm/mat.hpp>
//#include <boost/qvm/mat_traits.hpp>
//#include <boost/qvm/mat_access.hpp>
//#include <boost/qvm/mat_operations.hpp>

#include <boost/multiprecision/mpfi.hpp>

#include "../../include/ANTL/Cubic/generalFunctions.hpp"
#include "../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../include/ANTL/Cubic/CubicOrder.hpp"
#include "../../include/ANTL/Cubic/CubicElement.hpp"
#include "../../include/ANTL/Cubic/CubicIdeal.hpp"

#include "../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp"
#include "../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"

#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
using namespace NTL;
using namespace ANTL;
using NTL::ZZ;
using NTL::RR;
using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using boost::multiprecision::mpfi_float;
using std::cout;
using std::endl;
NTL_CLIENT

#endif
