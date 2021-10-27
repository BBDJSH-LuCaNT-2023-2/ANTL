#ifndef QUADRATICIDEALBASE_ZZ_TEST
#define QUADRATICIDEALBASE_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticIdealBase<ZZ>: Ideal Bases within different orders should not be equal", "[QuadraticIdealBase]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));
    QuadraticOrder<ZZ> quad_order2 = QuadraticOrder<ZZ>(ZZ(17));

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order2);

    REQUIRE(!quad_ideal_base1.IsEqual(quad_ideal_base2));

}

TEST_CASE("QuadraticIdealBase<ZZ>: Testing assign_prime()", "[QuadraticIdealBase]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));
    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign_prime(ZZ(17));

    REQUIRE(quad_ideal_base1.get_a() == ZZ(17));
    REQUIRE(quad_ideal_base1.get_b() == ZZ(9));
    REQUIRE(quad_ideal_base1.get_c() == ZZ(1));

}

#endif
