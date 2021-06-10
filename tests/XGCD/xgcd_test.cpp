#include <iostream>

// qvm matrix headers
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include "../../include/ANTL/XGCD/xgcd_plain.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>

using namespace NTL;

using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using std::cout;
using std::endl;
NTL_CLIENT

int main(){
long a,b,r,s,g;
a = 3;
b= 5;
XGCD(g, r,s,a,b);
cout << g << " = " << a << "*" << r <<  " + " << b << "*" << s << endl;


return 0;
}
