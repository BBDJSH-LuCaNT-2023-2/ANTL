#ifndef REDUCEPLAINREAL_ZZ_TEST
#define REDUCEPLAINREAL_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("ReducePlainReal<ZZ>: Ideals should be reduced upon calling reduce", "[ReducePlainReal]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(29));
    QuadraticNumber<ZZ> quad_number1 = QuadraticNumber<ZZ>(quad_order1);

    ReducePlainReal<ZZ> reduce_plain_real_object = ReducePlainReal<ZZ>();
    reduce_plain_real_object.set_RelativeGenerator(quad_number1);
    quad_order1.set_red_plain_real(reduce_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(5), ZZ(3), ZZ(-1));

    quad_order1.get_red_plain_real()->reduce(quad_ideal_base1);

    REQUIRE(quad_ideal_base1.is_reduced());

}

#endif
