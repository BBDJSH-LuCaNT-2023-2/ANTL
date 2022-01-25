#ifndef SQUARE_NUDUPL_ZZ_TEST
#define SQUARE_NUDUPL_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Square/SquareNudupl.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("SquareNudupl<ZZ>: Squaring an ideal", "[SquareNudupl]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(193));
    QuadraticNumber<ZZ> quad_number1 = QuadraticNumber<ZZ>(quad_order1);
    QuadraticNumber<ZZ> quad_number2 = QuadraticNumber<ZZ>(quad_order1);

    SquareNudupl<ZZ> sqr_nudupl_object = SquareNudupl<ZZ>();
    sqr_nudupl_object.set_RelativeGenerator(quad_number1);
    quad_order1.set_sqr_nudupl(sqr_nudupl_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    red_plain_real_object.set_RelativeGenerator(quad_number2);
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    QuadraticIdealBase<ZZ> quad_ideal_base_square = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(2), ZZ(13), ZZ(-3));

    sqr(quad_ideal_base_square, quad_ideal_base1);

    REQUIRE(quad_ideal_base_square.get_a() == 4);
    REQUIRE(quad_ideal_base_square.get_b() == 9);
    REQUIRE(quad_ideal_base_square.get_c() == -7);

    sqr(quad_ideal_base_square, quad_ideal_base_square);

    REQUIRE(quad_ideal_base_square.get_a() == 3);
    REQUIRE(quad_ideal_base_square.get_b() == 11);
    REQUIRE(quad_ideal_base_square.get_c() == -6);

}
#endif
