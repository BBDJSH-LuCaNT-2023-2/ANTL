#ifndef QUADRATICORDER_long_TEST
#define QUADRATICORDER_long_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticOrder.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticOrder<long>: Orders are equal iff discriminants are equal", "[QuadraticOrder]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));
    QuadraticOrder<long> quad_order2 = QuadraticOrder<long>(long(13));
    QuadraticOrder<long> quad_order3 = QuadraticOrder<long>(long(17));

    REQUIRE(quad_order1 == quad_order2);

    REQUIRE(quad_order1 != quad_order3);
}

TEST_CASE("QuadraticOrder<long>: get_discriminant returns the discriminant", "[QuadraticOrder]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));

    long expected_discriminant = long(13);

    REQUIRE(quad_order1.get_discriminant() == expected_discriminant);
}

TEST_CASE("QuadraticOrder<long>: Orders are real iff their discriminant is positive", "[QuadraticOrder]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));

    REQUIRE(quad_order1.is_real());

    REQUIRE_FALSE(quad_order1.is_imaginary());
}

#endif
