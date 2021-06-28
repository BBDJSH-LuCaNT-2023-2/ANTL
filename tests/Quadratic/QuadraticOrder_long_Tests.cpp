#ifndef QUADRATICORDER_ZZ_TEST
#define QUADRATICORDER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticOrder.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("Miscellaneous tests for QuadraticOrder<long>", "[capturing]") {

    long D1 = 13;
    long D2 = 17;

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(D1);
    QuadraticOrder<long> quad_order2 = QuadraticOrder<long>(D1);

    QuadraticOrder<long> quad_order3 = QuadraticOrder<long>(D2);

    long expected_discriminant = 13;

    REQUIRE(quad_order1 == quad_order2);

    REQUIRE(quad_order1 != quad_order3);

    REQUIRE(quad_order1.getDiscriminant() == expected_discriminant);

    REQUIRE(quad_order1.IsReal());

    REQUIRE_FALSE(quad_order1.IsImaginary());
}


#endif
