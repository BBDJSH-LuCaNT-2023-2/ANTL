#ifndef REDUCEFAST_ZZ_TEST
#define REDUCEFAST_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticIdealBase<ZZ>: Ideal Bases within different orders should not be equal", "[QuadraticIdealBase]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    ReduceFast<ZZ> reduce_fast_object = ReduceFast<ZZ>();
    quad_order1.set_red_fast(reduce_fast_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(9), ZZ(1));

    quad_order1.get_red_fast()->reduce(quad_ideal_base1);

    REQUIRE(quad_ideal_base1.IsReduced());

}

#endif


