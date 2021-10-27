#ifndef SQUARE_PLAIN_ZZ_TEST
#define SQUARE_PLAIN_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Square/SquarePlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("SquarePlain<ZZ>: Squaring two ideals", "[SquarePlain]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    SquarePlain<ZZ> sqr_plain_object = SquarePlain<ZZ>();
    quad_order1.set_sqr_plain(sqr_plain_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base_square = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(-9), ZZ(1));

    quad_order1.get_sqr_plain()->square(quad_ideal_base_square, quad_ideal_base1);

    REQUIRE(quad_ideal_base_square.get_a() == ZZ(289));
    REQUIRE(quad_ideal_base_square.get_b() == ZZ(59) );
    REQUIRE(quad_ideal_base_square.get_c() == ZZ(3)  );

}
#endif

