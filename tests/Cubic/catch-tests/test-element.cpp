#ifndef GUARD_test_element_cpp
#define GUARD_test_element_cpp

#include <cmath>
#include <ctime>
#include <boost/math/tools/polynomial.hpp>
#include "../../../include/ANTL/Cubic/CubicOrderNF.hpp"
#include "../../../include/ANTL/Cubic/CubicElementNF.hpp"
#include "../../../include/ANTL/Arithmetic/QQ.hpp"
// For floating point arithmetic error tolerance

const double DOUBLE_TOLERANCE = 0.0000001;


//using namespace NTL;
// 020-TestCase.cpp
// Let Catch provide main():
#include "../catch.hpp"




long coffee[3] = {1,2,3};

boost::math::tools::polynomial<long> poly1{{-1,-3,1,1}};
CubicOrderNF<long, double> Rufio(poly1);


TEST_CASE("CubicElement getter functions"){
  CubicElementNF<long, double> Ellie (&Rufio, coffee, 7);
  CubicElementNF<long, double> Est (&Rufio, coffee, 7);
  // This will be the constructor, the getters, and the assign function
  REQUIRE(Est.get_denom() == 7);
  REQUIRE(Est.get_u() == 1);
  REQUIRE(Est.get_x() == 2);
  REQUIRE(Est.get_y() == 3);
  REQUIRE(Est.is_equal(Ellie));

  SECTION("Reassignment and is zero"){
    long tea[3] = {0,0,0};
    CubicElementNF<long, double> Elias (&Rufio, tea, 1);
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

  coffee[0] = std::rand() % 1000; coffee[1] = 0; coffee[2] = 0;
  long bottom = (std::rand() % 1000) +1;

  SECTION("Trace of a rational number"){


    // random rational number
    QQ<long> rn1(coffee[0], denny);
    mul(rn1, rn1, 3);

    CubicElementNF<long, double> Est (&Rufio, coffee, 1);
    Est.trace(tempQ);
    REQUIRE(tempQ.isEqual((3*coffee[0])) );

    CubicElementNF<long, double> Ellie (&Rufio, coffee, denny);
    Ellie.trace(tempQ);
    REQUIRE(rn1.isEqual(tempQ));
  }

  SECTION("Trace of an arbitrary cubic number"){
    coffee[1] = std::rand() % 1000;
    coffee[2] = std::rand() % 1000;
    CubicElementNF<long, double>Ellie (&Rufio, coffee, bottom);

    Ellie.trace(tempQ);

    QQ<long> rn1( 3*coffee[0] - poly1[2]*coffee[1] -2*poly1[1]*coffee[2], bottom);
    REQUIRE(rn1.isEqual(tempQ));
  }
}

TEST_CASE("Addition"){
  CubicElementNF<long, double> Ellie (&Rufio, 2,2,2, 3);
  CubicElementNF<long, double> Est (&Rufio, 5, 1,112, 2);
  CubicElementNF<long, double> Erika (&Rufio, 0,0, 0, 5);
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
  CubicElementNF<long, double> Ellie (&Rufio, 2,2,2, 3);
  CubicElementNF<long, double> Est (&Rufio, 5, 1,112, 2);
  CubicElementNF<long, double> Erika (&Rufio, 0,0, 0, 5);
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
  CubicElementNF<long, double> Ellie (&Rufio, 2,2,2, 3);
  CubicElementNF<long, double> Est (&Rufio, 1, 0,0, 5);
  CubicElementNF<long, double> Erika (&Rufio, 0,0, 0, 5);
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

#endif
