#ifndef QUADRATICORDER_ZZ_TEST
#define QUADRATICORDER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticOrder.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticOrder<zz>: Orders are equal iff discriminants are equal", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));
    QuadraticOrder<ZZ> quad_order2 = QuadraticOrder<ZZ>(ZZ(13));
    QuadraticOrder<ZZ> quad_order3 = QuadraticOrder<ZZ>(ZZ(17));

    REQUIRE(quad_order1 == quad_order2);

    REQUIRE(quad_order1 != quad_order3);
}

TEST_CASE("QuadraticOrder<zz>: getDiscriminant returns the discriminant", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    ZZ expected_discriminant = ZZ(13);

    REQUIRE(quad_order1.getDiscriminant() == expected_discriminant);
}

TEST_CASE("QuadraticOrder<zz>: Orders are real iff their discriminant is positive", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    REQUIRE(quad_order1.IsReal());

    REQUIRE_FALSE(quad_order1.IsImaginary());
}

#endif
