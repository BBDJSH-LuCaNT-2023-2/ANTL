#ifndef SQUARE_NUDUPL_ZZ_TEST
#define SQUARE_NUDUPL_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Square/SquareNudupl.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("SquareNudupl<ZZ>: Squaring an ideal", "[SquareNudupl]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    SquareNudupl<ZZ> sqr_nudupl_object = SquareNudupl<ZZ>();
    quad_order1.set_sqr_nudupl(sqr_nudupl_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    QuadraticIdealBase<ZZ> quad_ideal_base_square = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(-9), ZZ(1));

    sqr(quad_ideal_base_square, quad_ideal_base1);

    REQUIRE(quad_ideal_base_square.get_a() == ZZ(1));
    REQUIRE(quad_ideal_base_square.get_b() == ZZ(3));
    REQUIRE(quad_ideal_base_square.get_c() == ZZ(-1));

}
#endif
