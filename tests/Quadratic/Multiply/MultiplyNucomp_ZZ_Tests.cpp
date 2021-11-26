#ifndef MULTIPLY_NUCOMP_ZZ_TEST
#define MULTIPLY_NUCOMP_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("MultiplyNucomp<ZZ>: Multiplying two ideals", "[MultiplyNucomp]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    MultiplyNucomp<ZZ> mul_nucomp_object = MultiplyNucomp<ZZ>();
    quad_order1.set_mul_nucomp(mul_nucomp_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order1);

    QuadraticIdealBase<ZZ> qib_pr = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base_product = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(-9), ZZ(1));
    quad_ideal_base2.assign(ZZ(17), ZZ(-9), ZZ(1));

    mul(quad_ideal_base_product, quad_ideal_base1, quad_ideal_base2);

    REQUIRE(quad_ideal_base_product.get_a() == ZZ(1));
    REQUIRE(quad_ideal_base_product.get_b() == ZZ(3));
    REQUIRE(quad_ideal_base_product.get_c() == ZZ(-1));

}
#endif

