#ifndef GUARD_test_order_cpp
#define GUARD_test_order_cpp

#include <cmath>
#include <ctime>
#include <boost/math/tools/polynomial.hpp>
#include "../../../include/ANTL/Cubic/CubicOrder.hpp"
#include "../../../include/ANTL/Cubic/CubicElement.hpp"
#include "../../../include/ANTL/Arithmetic/QQ.hpp"

#include "../catch.hpp"

// constants
//precision level requirement for double arithmetic
const double DOUBLE_TOLERANCE = 0.0000001;




boost::math::tools::polynomial<long> poly1{{4,1,1,1}};
CubicOrderComplex<long, double> Rufio(poly1);


TEST_CASE("Double: Cubic Order accessor functions"){
  REQUIRE(Rufio.get_discriminant() == -379);
  REQUIRE(Rufio.get_root1() - (-1.74295) < DOUBLE_TOLERANCE);
  REQUIRE(Rufio.get_rho2() - (1.29495) < DOUBLE_TOLERANCE);
  REQUIRE(Rufio.get_rho1() - (-1.74295) < DOUBLE_TOLERANCE);

  SECTION("equal tests"){
    boost::math::tools::polynomial<long> poly2{{4,1,1,1}};
    boost::math::tools::polynomial<long> poly3{{1,3,2,1}};
    CubicOrderComplex<long, double> ord2(poly2);
    CubicOrderReal<long, double> ord3(poly3);

    REQUIRE(Rufio.is_equal(ord2)==true);
    REQUIRE(ord2.is_equal(Rufio) == true);
    REQUIRE(Rufio.is_equal(ord3) == false);
    REQUIRE(ord3.is_equal(Rufio) == false);

    REQUIRE(is_equal(ord2, Rufio));
    REQUIRE(!is_equal(ord3, Rufio));
  }

  SECTION("root swapping"){
    boost::math::tools::polynomial<long> poly4{{-1,-3,1,1}};
    CubicOrderReal<long, double> ord4(poly4);

    REQUIRE(ord4.is_real());
    REQUIRE(ord4.get_root1() - (-2.170086) < DOUBLE_TOLERANCE);
    ord4.roots_swap_position(0,1);
    REQUIRE(ord4.get_root1() - (1.481194304) < DOUBLE_TOLERANCE);
    ord4.roots_swap_position(0,2);
    REQUIRE(ord4.get_root1() - (-0.3111078) < DOUBLE_TOLERANCE);
    //ord4.set_unit_strategy("BSGS");
    REQUIRE(ord4.get_fundamental_unit(0)->get_u() != 0);

  }

  SECTION("Ideal Multiplying"){
    REQUIRE(Rufio.get_fundamental_unit(0)->get_u() != 0);

  }
}

TEST_CASE("CubicOrderComplex functions"){

  SECTION("Test close_minimum"){
    boost::math::tools::polynomial<long> poly1{{11,5,5,1}};
    CubicOrder<ZZ, RR> * co_pointer = CubicOrder<ZZ, RR>::make_order(poly1);

    CubicIdeal<ZZ, RR> reduced_ideal1 = CubicIdeal<ZZ, RR>(NULL);
    CubicElement<ZZ, RR> minima       = CubicElement<ZZ, RR>(NULL);
    std::vector<RR> v1 = {RR(3.0)};
    CubicIdeal<ZZ,RR> ideal1 = CubicIdeal(co_pointer);
    std::cout << "just before" << ideal1.toString() << endl;
    co_pointer->close_minimum(reduced_ideal1,minima, ideal1, v1);

    // note that these results were verified using the COLLECT algorithm written in pari/gp
    REQUIRE(minima.get_u() == ZZ(3));
    REQUIRE(minima.get_x() == ZZ(-4));
    REQUIRE(minima.get_y() == ZZ(1));
    REQUIRE(minima.get_denom() == ZZ(1));
    REQUIRE(minima.get_order() == co_pointer);

    REQUIRE(reduced_ideal1.get_coeff(0,0) == ZZ(11));
    REQUIRE(reduced_ideal1.get_coeff(1,0) == ZZ(0));
    REQUIRE(reduced_ideal1.get_coeff(2,0) == ZZ(0));
    REQUIRE(reduced_ideal1.get_coeff(0,1) == ZZ(0));
    REQUIRE(reduced_ideal1.get_coeff(1,1) == ZZ(-6));
    REQUIRE(reduced_ideal1.get_coeff(2,1) == ZZ(1));
    REQUIRE(reduced_ideal1.get_coeff(0,2) == ZZ(11));
    REQUIRE(reduced_ideal1.get_coeff(1,2) == ZZ(-1));
    REQUIRE(reduced_ideal1.get_coeff(2,2) == ZZ(2));
    REQUIRE(reduced_ideal1.get_coeff(0,0) == ZZ(11));
    REQUIRE(reduced_ideal1.get_coeff(1,0) == ZZ(0));
    REQUIRE(reduced_ideal1.get_denom() == ZZ(1));
  }

}



#endif
