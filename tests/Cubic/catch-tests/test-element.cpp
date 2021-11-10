#ifndef GUARD_test_element_cpp
#define GUARD_test_element_cpp

#include <cmath>
#include <ctime>
#include <iostream>
#include <boost/math/tools/polynomial.hpp>
#include "../../../include/ANTL/Cubic/CubicOrder.hpp"
#include "../../../include/ANTL/Cubic/CubicElement.hpp"
#include "../../../include/ANTL/Arithmetic/QQ.hpp"
// For floating point arithmetic error tolerance

#include "../../../include/ANTL/Cubic/generalFunctions.hpp"
#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/CubicIdeal.hpp"

#include "../../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp"
#include "../../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"
const double DOUBLE_TOLERANCE = 0.0000001;

using NTL::ZZ;
using NTL::RR;
//using namespace NTL;
// 020-TestCase.cpp
// Let Catch provide main():
#include "../catch.hpp"




long elt_coeff[3] = {1,2,3};
boost::math::tools::polynomial<long> poly1{{-1,-3,1,1}};
CubicOrder<long, double> * Rufio = CubicOrder<long, double>::make_order(poly1);

ZZ newpoly[4];
ZZ elt_vector[3];


CubicOrder<ZZ, RR> * co_point;
CubicOrder<ZZ, RR> * second_order;


TEST_CASE("Constructors and ZZ,RR type tests"){
  newpoly[3] = 1;
  newpoly[2] = 1;
  newpoly[1] = -3;
  newpoly[0] = -1;
  polynomial<ZZ> const test_poly{{newpoly[0],newpoly[1],newpoly[2],newpoly[3] }};
  co_point = CubicOrder<ZZ, RR>::make_order(test_poly);
  CubicOrder<ZZ, RR> * Odie = co_point;
  ZZ vec2[4]; vec2[0] = 1; vec2[1] = 1; vec2[2] = 3; vec2[3] = 1;
  polynomial<ZZ> const poly2{{vec2[0],vec2[1],vec2[2],vec2[3] }};



  second_order = CubicOrder<ZZ, RR>::make_order(poly2);


  elt_vector[2] = 0;
  elt_vector[1] = 0;
  elt_vector[0] = 1;
  CubicElement<ZZ, RR> default_elt(Odie);
  CubicElement<ZZ, RR> elt1 (Odie, elt_vector, ZZ(1));
  CubicElement<ZZ, RR> elt2 (Odie, elt_vector[0], elt_vector[1],elt_vector[2], ZZ(7));

  SECTION("CubicElement getters, RR,ZZ"){
    REQUIRE(elt2.get_order() == Odie);

    REQUIRE(elt2.get_denom() == 7);
    REQUIRE(elt2.get_u() == 1);
    REQUIRE(elt2.get_x() == 0);
    REQUIRE(elt2.get_y() == 0);
    REQUIRE(elt1.is_equal(default_elt));
  }
      REQUIRE(elt2.get_order() == Odie);
  SECTION("CubicElement getters, RR,ZZ"){

    elt2.set_order(second_order);
    REQUIRE(elt2.get_order() == second_order);

    elt2.assign(elt1);
    elt2.set_order(elt1.get_order());
    REQUIRE(elt2.is_equal(elt1));
    elt2.assign(ZZ(0),ZZ(0),ZZ(0),ZZ(1));
    REQUIRE(elt2.is_zero() == true);
  }

}

TEST_CASE("Trace Function for ZZ and RR"){
  newpoly[3] = 1;
  newpoly[2] = 1;
  newpoly[1] = 4;
  newpoly[0] = 2;
  polynomial<ZZ> const test_poly{{newpoly[0],newpoly[1],newpoly[2],newpoly[3] }};
  co_point = CubicOrder<ZZ, RR>::make_order(test_poly);
  CubicOrder<ZZ, RR> * Odie = co_point;

  ZZ denny = to_ZZ(std::rand() % 1000);
  QQ<ZZ> tempQ(ZZ(1), ZZ(1));                     // placeholder rational
  ZZ eltZZ[3];
  eltZZ[0] = to_ZZ(std::rand() % 1000); eltZZ[1] = to_ZZ(0); eltZZ[2] = to_ZZ(0);
  ZZ bottom = ZZ((std::rand() % 1000) +1);

  SECTION("Trace of a rational number"){

    ZZ intermediate;
    // random rational number
    QQ<ZZ> rn1(eltZZ[0], denny);
    mul(rn1, rn1, ZZ(3));

    CubicElement<ZZ, RR> elt1 (Odie, eltZZ, ZZ(1));
    elt1.trace(tempQ);
    mul(intermediate, eltZZ[0], to<ZZ>(3));
    REQUIRE( tempQ.isEqual(intermediate));

    CubicElement<ZZ, RR> elt2 (Odie, eltZZ, denny);
    elt2.trace(tempQ);
    REQUIRE(rn1.isEqual(tempQ));

  }

  SECTION("Trace of an arbitrary cubic number"){
    eltZZ[1] = to<ZZ>(std::rand() % 1000);
    eltZZ[2] = to<ZZ>(std::rand() % 1000);
    CubicElement<ZZ, RR>elt2 (Odie, eltZZ, bottom);

    elt2.trace(tempQ);

    QQ<ZZ> rn1( ZZ(3)*eltZZ[0] - newpoly[2]*eltZZ[1] -ZZ(2)*newpoly[1]*eltZZ[2], bottom);
    REQUIRE(rn1.isEqual(tempQ));
  }

  SECTION("real value of an arbitrary cubic number"){
    eltZZ[0] = to<ZZ>(2);
    eltZZ[1] = to<ZZ>(16);
    eltZZ[2] = to<ZZ>(5);
    bottom = ZZ(1);
    CubicElement<ZZ, RR>elt2 (Odie, eltZZ, bottom);

    RR realholder;
    elt2.get_real_value(realholder);

    REQUIRE(realholder - (-7.77532584701453) < DOUBLE_TOL);
  }

  SECTION("real value of an arbitrary cubic number"){
    eltZZ[0] = to<ZZ>(2);
    eltZZ[1] = to<ZZ>(16);
    eltZZ[2] = to<ZZ>(5);
    bottom = ZZ(1);
    CubicElement<ZZ, RR> elt1 (Odie);
    CubicElement<ZZ, RR>elt2 (Odie, eltZZ, bottom);

    elt2.negate(elt1);

    REQUIRE(elt1.get_u() == to<ZZ>(-2));
    REQUIRE(elt1.get_x() == to<ZZ>(-16));
    REQUIRE(elt1.get_y() == to<ZZ>(-5));
  }
}














TEST_CASE("CubicElement getter functions (long, double)"){
  CubicOrder<long, double> * Rufio = CubicOrder<long, double>::make_order(poly1);
  CubicElement<long, double> Ellie (Rufio, elt_coeff, 7);
  CubicElement<long, double> Est (Rufio, elt_coeff, 7);
  // This will be the constructor, the getters, and the assign function
  REQUIRE(Est.get_denom() == 7);
  REQUIRE(Est.get_u() == 1);
  REQUIRE(Est.get_x() == 2);
  REQUIRE(Est.get_y() == 3);
  REQUIRE(Est.is_equal(Ellie));

  SECTION("Reassignment and is zero"){
    long tea[3] = {0,0,0};
    CubicElement<long, double> Elias (Rufio, tea, 1);
    Ellie.assign(Elias);

    REQUIRE(Ellie.get_u() == 0);
    REQUIRE(Ellie.get_x() == 0);
    REQUIRE(Ellie.get_y() == 0);
    REQUIRE(Ellie.is_zero() == true);

    Ellie.assign(3,3,3,1);

    REQUIRE(Ellie.get_u() == 3);
    REQUIRE(Ellie.get_x() == 3);
    REQUIRE(Ellie.get_y() == 3);
    REQUIRE(Ellie.get_denom() == 1);

    std::srand( (unsigned)std::time(NULL) );
    long rint = std::rand()%1000;
    Ellie.assign(rint);
    REQUIRE(Ellie.get_u() == rint);
    REQUIRE(Ellie.get_x() == 0);
    REQUIRE(Ellie.get_y() == 0);
    REQUIRE(Ellie.get_denom() == 1);

    QQ<long> gregory(6,12);
    Ellie.assign(gregory);
    REQUIRE(Ellie.get_u() == 1);
    REQUIRE(Ellie.get_x() == 0);
    REQUIRE(Ellie.get_y() == 0);
    REQUIRE(Ellie.get_denom() == 2);

  }

}

TEST_CASE("Trace Function"){
  long denny = std::rand() % 1000;
  QQ<long> tempQ(1, 1);                     // placeholder rational

  elt_coeff[0] = std::rand() % 1000; elt_coeff[1] = 0; elt_coeff[2] = 0;
  long bottom = (std::rand() % 1000) +1;

  SECTION("Trace of a rational number"){


    // random rational number
    QQ<long> rn1(elt_coeff[0], denny);
    mul(rn1, rn1, 3);

    CubicElement<long, double> Est (Rufio, elt_coeff, 1);
    Est.trace(tempQ);
    REQUIRE(tempQ.isEqual((3*elt_coeff[0])) );

    CubicElement<long, double> Ellie (Rufio, elt_coeff, denny);
    Ellie.trace(tempQ);
    REQUIRE(rn1.isEqual(tempQ));
  }

  SECTION("Trace of an arbitrary cubic number"){
    elt_coeff[1] = std::rand() % 1000;
    elt_coeff[2] = std::rand() % 1000;
    CubicElement<long, double>Ellie (Rufio, elt_coeff, bottom);

    Ellie.trace(tempQ);

    QQ<long> rn1( 3*elt_coeff[0] - poly1[2]*elt_coeff[1] -2*poly1[1]*elt_coeff[2], bottom);
    REQUIRE(rn1.isEqual(tempQ));
  }
}

TEST_CASE("Addition"){
  CubicElement<long, double> Ellie (Rufio, 2,2,2, 3);
  CubicElement<long, double> Est (Rufio, 5, 1,112, 2);
  CubicElement<long, double> Erika (Rufio, 0,0, 0, 5);
  QQ<long> tempQ(0, 2);
  SECTION("Add a constant"){
    add(Erika,Ellie, 0L);
    REQUIRE(Erika.is_equal(Ellie));

    add(Erika, Ellie, tempQ);
    REQUIRE(Erika.is_equal(Ellie));

    add(Erika, Ellie, 1005L);
    REQUIRE(Erika.get_u() == 3017);
    REQUIRE(Erika.get_x() == 2);
    REQUIRE(Erika.get_y() == 2);
    REQUIRE(Erika.get_denom() == 3);

    tempQ.setNumerator(7);
    tempQ.setDenominator(2);
    add(Erika, Ellie, tempQ);


    REQUIRE(Erika.get_u() == 25);
    REQUIRE(Erika.get_x() == 4);
    REQUIRE(Erika.get_y() == 4);
    REQUIRE(Erika.get_denom() == 6);

  }
  SECTION("Add two cubic elements"){
    add(Erika,Ellie, Est);

    REQUIRE(Erika.get_u() == 19);
    REQUIRE(Erika.get_x() == 7);
    REQUIRE(Erika.get_y() == 340);
    REQUIRE(Erika.get_denom() == 6);


  }

}

TEST_CASE("Subtraction"){
  CubicElement<long, double> Ellie (Rufio, 2,2,2, 3);
  CubicElement<long, double> Est (Rufio, 5, 1,112, 2);
  CubicElement<long, double> Erika (Rufio, 0,0, 0, 5);
  QQ<long> tempQ(0, 2);
  SECTION("Subtract a constant"){
    sub(Erika,Ellie, 0L);
    REQUIRE(Erika.is_equal(Ellie));

    sub(Erika, Ellie, tempQ);
    REQUIRE(Erika.is_equal(Ellie));

    sub(Erika, Ellie, 1005L);
    REQUIRE(Erika.get_u() == -3013);
    REQUIRE(Erika.get_x() == 2);
    REQUIRE(Erika.get_y() == 2);
    REQUIRE(Erika.get_denom() == 3);

    tempQ.setNumerator(7);
    tempQ.setDenominator(2);
    sub(Erika, Ellie, tempQ);


    REQUIRE(Erika.get_u() == -17);
    REQUIRE(Erika.get_x() == 4);
    REQUIRE(Erika.get_y() == 4);
    REQUIRE(Erika.get_denom() == 6);

  }
  SECTION("Sub two cubic elements"){
    sub(Erika,Ellie, Est);

    REQUIRE(Erika.get_u() == -11);
    REQUIRE(Erika.get_x() == 1);
    REQUIRE(Erika.get_y() == -332);
    REQUIRE(Erika.get_denom() == 6);


  }

}

TEST_CASE("Multiplication"){
  CubicElement<long, double> Ellie (Rufio, 2,2,2, 3);
  CubicElement<long, double> Est (Rufio, 1, 0,0, 5);
  CubicElement<long, double> Erika (Rufio, 0,0, 0, 5);
  QQ<long> tempQ(0, 2);

  SECTION("Mul by a constant"){
    mul(Erika,Ellie, 0L);
    REQUIRE(Erika.is_zero());

    mul(Erika,Ellie, tempQ);
    REQUIRE(Erika.is_zero());


    mul(Erika,Ellie, 150L);
    REQUIRE(Erika.get_u() == 100);
    REQUIRE(Erika.get_x() == 100);
    REQUIRE(Erika.get_y() == 100);
    REQUIRE(Erika.get_denom() == 1);

    tempQ.setNumerator(17);
    tempQ.setDenominator(11);
    mul(Erika,Ellie, tempQ);
    REQUIRE(Erika.get_u() == 34);
    REQUIRE(Erika.get_x() == 34);
    REQUIRE(Erika.get_y() == 34);
    REQUIRE(Erika.get_denom() == 33);
  }

  SECTION("Multiply two cubic elements"){
    // Note that the defining IBCF is x^3 + x^2 -3x -1

    // constant
    mul(Erika,Ellie, Est);

    REQUIRE(Erika.get_u() == 2);
    REQUIRE(Erika.get_x() == 2);
    REQUIRE(Erika.get_y() == 2);
    REQUIRE(Erika.get_denom() == 15);

    // just rho1
    Est.assign(0,1,0,2);
    mul(Erika,Ellie, Est);

    REQUIRE(Erika.get_u() == 1);
    REQUIRE(Erika.get_x() == 3);
    REQUIRE(Erika.get_y() == 1);
    REQUIRE(Erika.get_denom() == 3);

    // just rho2
    Est.assign(0,0,1,5);
    mul(Erika,Ellie, Est);

    REQUIRE(Erika.get_u() == 4);
    REQUIRE(Erika.get_x() == 8);
    REQUIRE(Erika.get_y() == 8);
    REQUIRE(Erika.get_denom() == 15);

    // All coeffs non-zero
    Est.assign(5,1,9,5);
    mul(Erika,Ellie, Est);

    REQUIRE(Erika.get_u() == 48);
    REQUIRE(Erika.get_x() == 88);
    REQUIRE(Erika.get_y() == 84);
    REQUIRE(Erika.get_denom() == 15);
  }
}

TEST_CASE("Division"){
  newpoly[3] = 1;
  newpoly[2] = 1;
  newpoly[1] = -2;
  newpoly[0] = -1;
  polynomial<ZZ> const test_poly{{newpoly[0],newpoly[1],newpoly[2],newpoly[3] }};
  co_point = CubicOrder<ZZ, RR>::make_order(test_poly);
  CubicOrder<ZZ, RR> * Odie = co_point;

  CubicElement<ZZ, RR> Ellie (Odie, ZZ(5),ZZ(0),ZZ(0), ZZ(2));
  CubicElement<ZZ, RR> Est (Odie, ZZ(1), ZZ(0),ZZ(0), ZZ(1));
  CubicElement<ZZ, RR> Erika (Odie, ZZ(2),ZZ(0), ZZ(0), ZZ(5));
  QQ<ZZ> tempQ(ZZ(0), ZZ(2));

  SECTION("Inverse test"){
    mul(Erika, Ellie, Erika);
    std::cout << Erika.toString() << std::endl;
    REQUIRE(Erika.is_equal(Est));

    Ellie.inverse(Erika);
    std::cout << Ellie.toString() << std::endl;
    mul(Erika, Ellie, Erika);
    REQUIRE(Erika.is_equal(Est));

    Ellie.assign(ZZ(14), ZZ(3), ZZ(7), ZZ(2));
    Ellie.inverse(Erika);
    mul(Erika, Ellie, Erika);
    REQUIRE(Erika.is_equal(Est));

    div(Erika, Ellie, Ellie);
    REQUIRE(Erika.is_equal(Est));

    Erika.assign(ZZ(1), ZZ(1), ZZ(0), ZZ(1));
    REQUIRE(Erika.get_u() ==ZZ(1));
    REQUIRE(Erika.get_x() == ZZ(1));
    REQUIRE(Erika.get_y() == ZZ(0));
    REQUIRE(Erika.get_denom() == ZZ(1));
    Ellie.inverse(Erika);
    mul(Ellie, Ellie, Erika);
    //std::cout << Erika.toString() << std::endl;
    REQUIRE(Ellie.get_u() == Est.get_u());
    REQUIRE(Ellie.get_x() == Est.get_x());
    REQUIRE(Ellie.get_y() == Est.get_y());
    REQUIRE(Ellie.get_denom() == Est.get_denom());
    //REQUIRE(Erika.is_equal(Est));
  }


  //SECTION("Multiply two cubic elements"){
//
  //}
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////








#endif
