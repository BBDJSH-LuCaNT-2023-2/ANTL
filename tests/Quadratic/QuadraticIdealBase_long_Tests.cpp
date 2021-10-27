#ifndef QUADRATICIDEALBASE_long_TEST
#define QUADRATICIDEALBASE_long_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticIdealBase<long>: Ideal Bases within different orders should not be equal", "[QuadraticIdealBase]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));
    QuadraticOrder<long> quad_order2 = QuadraticOrder<long>(long(17));

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base2 = QuadraticIdealBase<long>(quad_order2);

    REQUIRE(!quad_ideal_base1.IsEqual(quad_ideal_base2));

}

TEST_CASE("QuadraticIdealBase<long>: Testing assign_prime()", "[QuadraticIdealBase]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));
    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign_prime(17);

    REQUIRE(quad_ideal_base1.get_a() == 17);
    REQUIRE(quad_ideal_base1.get_b() == 9);
    REQUIRE(quad_ideal_base1.get_c() == 1);

}

#endif
