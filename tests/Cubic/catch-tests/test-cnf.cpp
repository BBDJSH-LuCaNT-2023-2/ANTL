#ifndef GUARD_testcnf_cpp
#define GUARD_testcnf_cpp

#include <cmath>
#include <boost/math/tools/polynomial.hpp>
#include "../../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../../include/ANTL/Arithmetic/QQ.hpp"
// For floating point arithmetic error tolerance
const double DOUBLE_TOLERANCE = 0.0000001;


using namespace NTL;
using namespace ANTL;
using boost::math::tools::polynomial;
// 010-TestCase.cpp

//using namespace NTL;

// Let Catch provide main():

#include "../catch.hpp"



  int cubic_type;
  long a,b,c,d;

// This test case is for GeneralTemplateFunctions
TEST_CASE("Root-finding function SolveP3, monic, double") {


  double rootArray[3];
  SECTION("Integral Roots"){
    a = 1; b = -2; c = -1; d = 2;
    polynomial<long> P1{{d,c,b,a}};


    cubic_type = cardano<long, double>(P1,rootArray);

    REQUIRE ( ANTL::abs(rootArray[0] - (-1.0) ) < DOUBLE_TOLERANCE );
    REQUIRE ( ANTL::abs(rootArray[2] - 1.0) < DOUBLE_TOLERANCE );
    REQUIRE ( ANTL::abs(rootArray[1] - 2.0) <DOUBLE_TOLERANCE );
  }

  SECTION("Complex Roots"){
    a = 1; b = 1; c = 1; d = 4;
    polynomial<long> P1{{d,c,b,a}};
    cubic_type = cardano<long, double>(P1,rootArray);

    REQUIRE ( ANTL::abs(rootArray[0] - (-1.7429592) ) < DOUBLE_TOLERANCE );
    REQUIRE ( ANTL::abs(rootArray[1] - (0.3714796) ) < DOUBLE_TOLERANCE );
    REQUIRE ( ANTL::abs(ANTL::abs(rootArray[2]) - (1.4686560) ) < DOUBLE_TOLERANCE );
  }

  SECTION("Not Monic, real"){
    a = 2; b = 1; c = -5; d = -2; polynomial<long> P1{{d,c,b,a}};
    cubic_type = cardano<long, double>(P1, rootArray);

    REQUIRE ( ANTL::abs(rootArray[0] - (1.5419361 ) ) < DOUBLE_TOLERANCE );
    REQUIRE ( ANTL::abs(rootArray[1] - (-1.6485352) ) < DOUBLE_TOLERANCE );
    REQUIRE ( ANTL::abs(rootArray[2] - (-0.3934009) ) < DOUBLE_TOLERANCE );
  }
} // end test case


boost::math::tools::polynomial<long> poly1{{-1,-3,1,1}};
boost::math::tools::polynomial<boost::multiprecision::mpz_int> poly2{{-1,-3,1,1}};
TEST_CASE("Discriminant function"){


  SECTION("Binary Cubic Form, discriminant long and mpz_int"){
    REQUIRE(discriminant_bcf(poly1) == 148);
    REQUIRE(discriminant_bcf(poly2) == 148);
  }
}


/*
// Here are tests for the RealCubicNF
TEST_CASE("Constructor of RealCubicNF, and accessors"){
  RealCubicNumberField<long, double> Rufio(poly1);
  boost::math::tools::polynomial<long> * Rpoly;
  Rpoly = Rufio.get_defining_polynomial();

  SECTION("Testing Global accessors and RealCubicNF specific methods"){
    REQUIRE(Rufio.get_discriminant() == 148);
    REQUIRE( Rpoly->degree() == 3 );
    REQUIRE( (*Rpoly)[0] == -1 );
    REQUIRE( (*Rpoly)[1] == -3 );
    REQUIRE( (*Rpoly)[2] == 1 );
    REQUIRE(Rufio.is_real() == true);
    REQUIRE(Rufio.is_complex() == false);
    REQUIRE(ANTL::abs(Rufio.get_root(0)) > 1.0);

    double rtemp[3] = {Rufio.get_root(0),Rufio.get_root(1), Rufio.get_root(2)};
    Rufio.roots_swap_position(0,2);
    REQUIRE(rtemp[0] == Rufio.get_root(2));
    REQUIRE(rtemp[1] == Rufio.get_root(1));
    REQUIRE(rtemp[2] == Rufio.get_root(0));

  }
}*/

/*
TEST_CASE("Constructor for ComplexCubicNumberField and accessors"){

  poly1 = {3,1,1,1};
  ComplexCubicNumberField<long, double> Complejo(poly1);

  SECTION("Correctness of polynomial properties"){
    REQUIRE(poly1[0] == 3);
    REQUIRE(poly1[1] == 1);
    REQUIRE(poly1[2] == 1);
    REQUIRE(poly1[3] == 1);
    REQUIRE(Complejo.get_discriminant() == -204);
  }

  SECTION("ComplexCubicNumberField methods"){

    REQUIRE(Complejo.is_complex() == true);
    REQUIRE(Complejo.is_real() == false);
    REQUIRE(Complejo.get_real_root() - (-1.574743) < DOUBLE_TOLERANCE);
    REQUIRE(Complejo.get_complex_root().real() - (0.2873715) < DOUBLE_TOLERANCE);
    REQUIRE(ANTL::abs(Complejo.get_complex_root().imag()) - (1.3499963) < DOUBLE_TOLERANCE);
  }
}
*/


#endif
