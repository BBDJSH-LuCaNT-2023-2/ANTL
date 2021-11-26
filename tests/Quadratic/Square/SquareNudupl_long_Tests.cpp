#ifndef SQUARE_NUDUPL_LONG_TEST
#define SQUARE_NUDUPL_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Square/SquareNudupl.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("SquareNudupl<long>: Squaring an ideal", "[SquareNudupl]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));

    SquareNudupl<long> sqr_nudupl_object = SquareNudupl<long>();
    quad_order1.set_sqr_nudupl(sqr_nudupl_object);

    ReducePlainReal<long> red_plain_real_object = ReducePlainReal<long>();
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);

    QuadraticIdealBase<long> quad_ideal_base_square = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, -9, 1);

    sqr(quad_ideal_base_square, quad_ideal_base1);

    REQUIRE(quad_ideal_base_square.get_a() == 1);
    REQUIRE(quad_ideal_base_square.get_b() == 3);
    REQUIRE(quad_ideal_base_square.get_c() == -1);

}
#endif


