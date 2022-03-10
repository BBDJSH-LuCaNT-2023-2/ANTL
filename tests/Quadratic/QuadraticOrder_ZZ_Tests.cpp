#ifndef QUADRATICORDER_ZZ_TEST
#define QUADRATICORDER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticOrder.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticOrder<ZZ>: Orders are equal iff discriminants are equal", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1{ZZ(13)};
    QuadraticOrder<ZZ> quad_order2{ZZ(13)};
    QuadraticOrder<ZZ> quad_order3{ZZ(17)};

    REQUIRE(quad_order1 == quad_order2);

    REQUIRE(quad_order1 != quad_order3);
}

TEST_CASE("QuadraticOrder<ZZ>: getDiscriminant returns the discriminant", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1{ZZ(13)};

    ZZ expected_discriminant = ZZ(13);

    REQUIRE(quad_order1.getDiscriminant() == expected_discriminant);
}

TEST_CASE("QuadraticOrder<ZZ>: Orders are real iff their discriminant is positive", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1{ZZ(13)};

    REQUIRE(quad_order1.IsReal());

    REQUIRE_FALSE(quad_order1.IsImaginary());
}

TEST_CASE("QuadraticOrder<ZZ>: Correctness of Lfunc related functions", "[QuadraticOrder]") {

    QuadraticOrder<ZZ> quad_order1{ZZ(10121)};

    ZZ hREstimate = quad_order1.approximate_hR();
    //ZZ hRLowerBound = quad_order1.lower_bound_hR();
    ZZ hREstimateError = quad_order1.estimate_hR_error();

    std::cout << "hREstimate is " << hREstimate << std::endl;
    //REQUIRE(hRLowerBound < hREstimate);
    std::cout << "hREstimateError is " << hREstimateError << std::endl;
    REQUIRE(hREstimate < 5000);

}

#endif
