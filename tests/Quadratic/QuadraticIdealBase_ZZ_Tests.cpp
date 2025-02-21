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

TEST_CASE("QuadraticIdealBase<ZZ>: Testing is_reduced()", "[QuadraticIdealBase]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(5), ZZ(3), ZZ(-1)); //A form which is not reduced [BV07, pg. 109]
    quad_ideal_base2.assign(ZZ(1), ZZ(5), ZZ(-1)); //A form which is reduced     [BV07, pg. 109]

    REQUIRE(quad_ideal_base1.is_reduced() == false);
    REQUIRE(quad_ideal_base2.is_reduced() == true);
}

TEST_CASE("QuadraticIdealBase<ZZ>: Testing is_normal()", "[QuadraticIdealBase]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(5), ZZ(5), ZZ(1));  //A form which is normal [BV07, pg. 107]
    quad_ideal_base2.assign(ZZ(-3), ZZ(5), ZZ(4)); //A form which is normal [BV07, pg. 107]

    REQUIRE(quad_ideal_base1.is_normal() == true);
    REQUIRE(quad_ideal_base2.is_normal() == true);
}

#endif
