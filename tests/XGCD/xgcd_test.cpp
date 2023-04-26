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
        XGCD_PLAIN(g,x,y,a,b);
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

TEMPLATE_TEST_CASE("XGCD_LEFT_PLAIN tests","[XGCD][XGCD_LEFT][XGCD_LEFT_PLAIN]", int64_t, ZZ){
    TestType a, b, u, g, expected_g, expected_u;
    a = 0;
    b = 0;
//    XGCDLeftPlainTestInstance inst(a,b);

    SECTION("Basic Test"){ // (5,3) = 1, specifically -1*5 + 2*3 = 1
        a = 5;
        b = 3;
        expected_u = -1;
        expected_g = 1;
        cout << "Starting XGCD_LEFT_PLAIN test" << endl; 
        XGCD_LEFT_PLAIN(g,u,a,b);
        REQUIRE(g==expected_g);
        REQUIRE(u==expected_u);
    }

    SECTION("Negative Test"){ // (-6, -10) = 2, specifically -2*-6 + 1*-10 = 2
        a = -6;
        b = -10;
        expected_g=2;
        expected_u=-2;
        XGCD_LEFT_PLAIN(g,u,a,b);
        cout << g << " " << u << " " << a << " _ " << b << endl;
        cout << (g-a*u)/b << endl;
        REQUIRE(g==expected_g);
        //REQUIRE(u==expected_u);

    }
}

TEMPLATE_TEST_CASE("XGCD_BINARY_L2R tests", "[XGCD][XGCD_PLAIN][XGCD_BINARY_L2R]", int64_t){
    TestType g,x,y,b,a, expected_g;
    //vector<TestType> sol;
    //XGCDBinaryL2RPlainTestInstance<TestType>* inst = new XGCDBinaryL2RPlainTestInstance<TestType>(a,b);
    SECTION("Basic Test"){
        a = 3;
        b = 5;
        expected_g=1;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g==expected_g);
    }
    SECTION("Basic Test Reversed"){
        a = 5;
        b = 3;
        expected_g = 1;
        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g==expected_g);
    }
    SECTION("Common Factor Test"){
        //common factor test: (5,20) = 5
        a = 5;
        b = 20;
        expected_g = 5;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);
    }

    SECTION("Unit Input Test"){
        //unit test: (1,999) = 1
        a = 1;
        b = 999;
        expected_g = 1;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);
    }

    SECTION("Negative Input Tests"){
        //unit test: (+/-3, +/-6) = 3
        a = -3;
        b = 6;
        expected_g = 3;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);

        a = 3;
        b = -6;
        expected_g = 3;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);

        a = -3;
        b = 6;
        expected_g = 3;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);
    }

    SECTION("64-bit Random Positive Test"){
        //arbitrary large (2^32 < a,b < 2^63) positive test: (3166167471260038366, 2078992898117306689) = 1
        a = 3166167471260038366;
        b = 2078992898117306689;
        expected_g = 1;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);

    }

    SECTION("64-bit Random Negative Test"){
        //arbitrary large (-2^32 > a,b > 2^63) negative test: (-3867470587490682194, -6531477986582055176) = 2
        a = -3867470587490682194;
        b =  -6531477986582055176;
        expected_g = 2;

        XGCD_BINARY_L2R(g,x,y,a,b);
        REQUIRE(g == expected_g);

    }
    
}

TEMPLATE_TEST_CASE("XGCD_BINARY_L2R_LEFT", "[XGCD][XGCD_LEFT][XGCD_BINARY_L2R]", int64_t){
    TestType g,x=2,b,a, expected_g, expected_x;
    SECTION("Basic Tests"){
        a = 3;
        b = 5;
        expected_g = 1;
        expected_x = 2;
        XGCD_BINARY_L2R_LEFT(g, x, a, b);
        REQUIRE(g==expected_g);
        REQUIRE(x==expected_x);

        a = 5;
        b = 3;
        expected_g = 1;
        expected_x = -1;
        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(g==expected_g);
        REQUIRE(x==expected_x);
    }
    SECTION("Common Factor Test"){
        //common factor test: (5,20) = 5
        a = 5;
        b = 20;
        expected_g = 5;
        expected_x = 1;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(g == expected_g);
        REQUIRE(x == expected_x);
    }

    SECTION("Unit Input Test"){
        //unit test: (1,999) = 1
        a = 1;
        b = 999;
        expected_g = 1;
        expected_x = 1;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(g == expected_g);
        REQUIRE(x == expected_x);
    }

    SECTION("Negative Input Tests"){
        //unit test: (+/-3, +/-6) = 3
        a = -3;
        b = 6;
        expected_g = 3;
        expected_x = -1;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(x == expected_x);
        REQUIRE(g == expected_g);

        a = 3;
        b = -6;
        expected_g = 3;
        expected_x = 1;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(g == expected_g);
        REQUIRE(x == expected_x);

        a = -3;
        b = 6;
        expected_g = 3;
        expected_x = -1;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(x == expected_x);
        REQUIRE(g == expected_g);
    }

    SECTION("64-bit Random Positive Test"){
        //arbitrary large (2^32 < a,b < 2^63) positive test: (3166167471260038366, 2078992898117306689) = 1
        a = 3166167471260038366;
        b = 2078992898117306689;
        expected_g = 1;
        expected_x = -918826731003719737;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(g == expected_g);
        REQUIRE(x == expected_x);

    }

    SECTION("64-bit Random Negative Test"){
        //arbitrary large (-2^32 > a,b > 2^63) negative test: (-3867470587490682194, -6531477986582055176) = 2
        a = -3867470587490682194;
        b =  -6531477986582055176;
        expected_g = 2;
        expected_x = -357142171114937269;

        XGCD_BINARY_L2R_LEFT(g,x,a,b);
        REQUIRE(g == expected_g);

    }

}

TEMPLATE_TEST_CASE("XGCD_BINARY_L2R_PARTIAL", "[XGCD][XGCD_PARTIAL][XGCD_BINARY_L2R]", int64_t){
    //TODO: Add partial tests
    TestType z, R1, R2, C1, C2;
    TestType bound;
    SECTION("Basic Test"){
        //a=5850, b=3596
        //r2=R_{-1}=b
        //r1=R_0 = a-floor(a/b)b
        TestType a = 5850;
        TestType b = 3596;
        R2 = b;
        R1 = (int64_t)a-floor(a/b)*b;
        cout << "R1 precompute: " << R1 << endl;
        bound=50;
        XGCD_PARTIAL_BINARY_L2R(z,R2,R1, C2,C1,bound);
        cout << "Z: " << z << " R2 " << R2 << " R1 " << R1 << " C2 " << C2 << " C1 " << C1 << endl;
        REQUIRE(R1>=R2);
        REQUIRE(R1<=bound);
        REQUIRE(bound < R2);   
    }

}
