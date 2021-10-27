#ifndef REDUCEFAST_LONG_TEST
#define REDUCEFAST_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticIdealBase<long>: Ideal Bases within different orders should not be equal", "[QuadraticIdealBase]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(13);

    ReduceFast<long> reduce_fast_object = ReduceFast<long>();
    quad_order1.set_red_fast(reduce_fast_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, 9, 1);

    quad_order1.get_red_fast()->reduce(quad_ideal_base1);

    REQUIRE(quad_ideal_base1.IsReduced());

}

#endif


