#ifndef GUARD_test_order_cpp
#define GUARD_test_order_cpp

#include <cmath>
#include <ctime>
#include <boost/math/tools/polynomial.hpp>
#include "../../../include/ANTL/Cubic/CubicOrderNF.hpp"
#include "../../../include/ANTL/Cubic/CubicElementNF.hpp"
#include "../../../include/ANTL/Arithmetic/QQ.hpp"

#include "../catch.hpp"

// constants
//precision level requirement for double arithmetic
const double DOUBLE_TOLERANCE = 0.0000001;




boost::math::tools::polynomial<long> poly1{{4,1,1,1}};
CubicOrderNF<long, double> Rufio(poly1);


TEST_CASE("Double: Cubic Order accessor functions"){
  REQUIRE(Rufio.get_discriminant() == -379);
  REQUIRE(Rufio.get_root1() - (-1.74295) < DOUBLE_TOLERANCE);
  REQUIRE(Rufio.get_rho2() - (1.29495) < DOUBLE_TOLERANCE);
  REQUIRE(Rufio.get_rho1() - (-1.74295) < DOUBLE_TOLERANCE);

  SECTION("equal tests"){
    boost::math::tools::polynomial<long> poly2{{4,1,1,1}};
    boost::math::tools::polynomial<long> poly3{{1,3,2,1}};
    CubicOrderNF<long, double> ord2(poly2);
    CubicOrderNF<long, double> ord3(poly3);

    REQUIRE(Rufio.is_equal(ord2)==true);
    REQUIRE(ord2.is_equal(Rufio) == true);
    REQUIRE(Rufio.is_equal(ord3) == false);
    REQUIRE(ord3.is_equal(Rufio) == false);

    REQUIRE(is_equal(ord2, Rufio));
    REQUIRE(!is_equal(ord3, Rufio));
  }

  SECTION("root swapping"){
    boost::math::tools::polynomial<long> poly4{{-1,-3,1,1}};
    CubicOrderNF<long, double> ord4(poly4);

    REQUIRE(ord4.get_root1() - (-2.170086) < DOUBLE_TOLERANCE);
    ord4.roots_swap_position(1,2);
    REQUIRE(ord4.get_root1() - (-0.311108) < DOUBLE_TOLERANCE);
    ord4.roots_swap_position(1,3);
    REQUIRE(ord4.get_root1() - (1.48119430) < DOUBLE_TOLERANCE);

  }

  SECTION("Ideal Multiplying"){
    
  }
}



#endif
