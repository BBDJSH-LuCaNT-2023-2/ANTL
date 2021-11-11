#ifndef MULTIPLY_NUCOMP_LONG_TEST
#define MULTIPLY_NUCOMP_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("MultiplyNucomp<long>: Multiplying two ideals", "[MultiplyNucomp]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(13);

    MultiplyNucomp<long> mul_nucomp_object = MultiplyNucomp<long>();
    quad_order1.set_mul_nucomp(mul_nucomp_object);

    ReducePlainReal<long> red_plain_real_object = ReducePlainReal<long>();
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base2 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base_product = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, -9, 1);
    quad_ideal_base2.assign(17, -9, 1);

    mul(quad_ideal_base_product, quad_ideal_base1, quad_ideal_base2);

    REQUIRE(quad_ideal_base_product.get_a() == long(1));
    REQUIRE(quad_ideal_base_product.get_b() == long(3));
    REQUIRE(quad_ideal_base_product.get_c() == long(-1));

}
#endif


