#ifndef REDUCEPLAINREAL_LONG_TEST
#define REDUCEPLAINREAL_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("ReducePlainReal<long>: Ideals should be reduced upon calling reduce", "[ReducePlainReal]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(13);

    ReducePlainReal<long> reduce_plain_real_object = ReducePlainReal<long>();
    quad_order1.set_red_plain_real(reduce_plain_real_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, 9, 1);

    quad_order1.get_red_plain_real()->reduce(quad_ideal_base1);

    REQUIRE(quad_ideal_base1.IsReduced());

}

#endif


