#ifndef SQUARE_PLAIN_LONG_TEST
#define SQUARE_PLAIN_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Square/SquarePlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("SquarePlain<long>: Squaring two ideals", "[SquarePlain]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(13);

    SquarePlain<long> sqr_plain_object = SquarePlain<long>();
    quad_order1.set_sqr_plain(sqr_plain_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base_square = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(long(17), long(-9), long(1));

    quad_order1.get_sqr_plain()->square(quad_ideal_base_square, quad_ideal_base1);

    REQUIRE(quad_ideal_base_square.get_a() == 289);
    REQUIRE(quad_ideal_base_square.get_b() == 59 );
    REQUIRE(quad_ideal_base_square.get_c() == 3  );

}
#endif

