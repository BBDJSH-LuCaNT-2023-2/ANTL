#include <iostream>

#include "tests/catch.hpp"

// qvm matrix headers
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include "../../include/ANTL/XGCD/xgcd_plain.hpp"
#include "../../include/ANTL/XGCD/xgcd_binary_l2r.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <vector>

using namespace NTL;

using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using std::cout;
using std::endl;
using std::vector;
NTL_CLIENT

TEMPLATE_TEST_CASE("XGCD_PLAIN tests", "[XGCD][XGCD_PLAIN]", int64_t){
    TestType a,b,x,y,g, expected_g;
    a = 0;
    b = 0;
    //XGCDPlainTestInstance<TestType>* inst = new XGCDPlainTestInstance<TestType>(a,b);

    SECTION("Basic Test"){
        //basic initial test: (3,5) = 1
        a = 3;
        b = 5;
        expected_g = 1;
        cout << "A " << a << " B " << b << " G " << g << endl;
        XGCD_PLAIN(g,x,y,a,b);
        cout << "A " << a << " B " << b << " G " << g << endl;
        REQUIRE(expected_g == g);
    }
    SECTION("Reverse Basic Test"){
        a = 5;
        b = 3;
        expected_g = 1;
    
        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(expected_g == g);
    }

    SECTION("Common Factor Test"){
        //common factor test: (5,20) = 5
        a = 5;
        b = 20;
        expected_g = 5;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);
    }

    SECTION("Unit Input Test"){
        //unit test: (1,999) = 1
        a = 1;
        b = 999;
        expected_g = 1;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);
    }

    SECTION("Negative Input Test"){
        //unit test: (+/-3, +/-6) = 3
        a = -3;
        b = 6;
        expected_g = 3;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);

        a = 3;
        b = -6;
        expected_g = 3;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);

        a = -3;
        b = 6;
        expected_g = 3;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);
    }

    SECTION("Large Random Positive Test"){
        //arbitrary large (2^32 < a,b < 2^63) positive test: (3166167471260038366, 2078992898117306689) = 1
        a = 3166167471260038366;
        b = 2078992898117306689;
        expected_g = 1;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);

    }

    SECTION("Large Random Negative Test"){
        //arbitrary large (-2^32 > a,b > 2^63) negative test: (-3867470587490682194, -6531477986582055176) = 2
        a = -3867470587490682194;
        b =  -6531477986582055176;
        expected_g = 2;

        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);

    }
    SECTION("0 Input Test"){
        //0 value tests
        a = 0;
        b = 100;
        expected_g = 100;
        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);

        a = 100;
        b = 0;
        expected_g = 100;
        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);
        
        a = 0;
        b = 0;
        expected_g = 0;
        XGCD_PLAIN(g,x,y,a,b);
        REQUIRE(g == expected_g);
        
    }



}

/*TEMPLATE_TEST_CASE("XGCD_LEFT_PLAIN tests","[XGCD][XGCD_LEFT][XGCD_LEFT_PLAIN]", int64_t){
    TestType a, b, u, g, expected_g, expected_u;
    a = 0;
    b = 0;
//    XGCDLeftPlainTestInstance inst(a,b);

    SECTION("Basic Test"){ // (3,5) = 1, specifically 2*3 -1*5 = 1
        a = 3;
        b = 5;
        expected_u = 2;
        expected_g = 1;

        XGCD_LEFT_PLAIN(g,u,a,b);
        REQUIRE(g==expected_g);
        REQUIRE(u==expected_u);

    }
    SECTION("Reverse Basic Test"){
        a = 5;
        b = 3;
        expected_u = -1;
        expected_g = 1;

        XGCD_LEFT_PLAIN(g,u,a,b);
        REQUIRE(g==expected_g);
        REQUIRE(u==expected_u);
    }
}
*/
TEMPLATE_TEST_CASE("XGCD_BINARY_L2R tests", "[XGCD][XGCD_BINARY_L2R]", int64_t){
    TestType g,x,y,b,a, expected_g;
    //vector<TestType> sol;
    //XGCDBinaryL2RPlainTestInstance<TestType>* inst = new XGCDBinaryL2RPlainTestInstance<TestType>(a,b);
    SECTION("Basic Test"){
        x = 1;
        y = 1;
        a = 3;
        b = 5;
        g = 2;
        expected_g=1;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g==expected_g);
    }
}