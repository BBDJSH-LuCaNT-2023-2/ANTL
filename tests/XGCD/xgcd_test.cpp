#include <iostream>

#include "tests/catch.hpp"

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

/*int main(){
long a,b,r,s,g;
a = 3;
b= 5;
XGCD(g, r,s,a,b);
cout << g << " = " << a << "*" << r <<  " + " << b << "*" << s << endl;

XGCD_PLAIN(r,s,a,b);
cout << r << " " << s << " " << a << " " << b << "\n";
return 0;
}*/

TEST_CASE("[XGCD]: int64_t PLAIN"){
    int64_t g,a,b,r,s;

    //basic initial test: (3,5) = 1
    a = 3;
    b = 5;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 1);

    //reverse magnitude test: (5,3) = 1
    a = 5;
    b = 3;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 1);


    //common factor test: (5,20) = 5
    a = 5;
    b = 20;
    XGCD_PLAIN(g, r,s,a,b);
    REQUIRE(g == 5);

    //unit test: (1, 999) = 1
    a = 1;
    b = 999;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 1);

    //sign permutation tests: (+/- 3, +/- 6)
    a = 3;
    b = 6;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 3);

    a = -3;
    b = 6;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 3);

    a = 3;
    b = -6;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 3);

    a = -3;
    b = -6;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 3);

    //arbitrary large (2^32 < a,b < 2^63) positive test: (3166167471260038366, 2078992898117306689) = 1
    a = 3166167471260038366;
    b = 2078992898117306689;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 1);

    //arbitrary large (-2^32 > a,b > 2^63) negative test: (-3867470587490682194, -6531477986582055176) = 2 
    a = -3867470587490682194;
    b = -6531477986582055176;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g == 2);

    a=79;
    b=420;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g==1);


    //0 value tests
    a = 0;
    b = 100;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g==100);

    a = 50;
    b = 0;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g==50);

    a = 0;
    b = 0;
    XGCD_PLAIN(g,r,s,a,b);
    REQUIRE(g==0);
    //cout g << " " << r << " " << s << " " << a << " " << b << "\n";




}
