#ifndef MULTIPLY_NUCOMP_ZZ_TEST
#define MULTIPLY_NUCOMP_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("MultiplyNucomp<ZZ>: Multiplying two ideals", "[MultiplyNucomp]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(193));
    QuadraticNumber<ZZ> quad_number1 = QuadraticNumber<ZZ>(quad_order1);
    QuadraticNumber<ZZ> quad_number2 = QuadraticNumber<ZZ>(quad_order1);

    MultiplyNucomp<ZZ> mul_nucomp_object = MultiplyNucomp<ZZ>();
    mul_nucomp_object.set_RelativeGenerator(quad_number1);
    quad_order1.set_mul_nucomp(mul_nucomp_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    red_plain_real_object.set_RelativeGenerator(quad_number2);
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order1);

    QuadraticIdealBase<ZZ> qib_pr = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base_product = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(2), ZZ(13), ZZ(-3));
    quad_ideal_base2.assign(ZZ(4), ZZ(9), ZZ(-7));

    mul(quad_ideal_base_product, quad_ideal_base1, quad_ideal_base2);

    //REQUIRE(quad_ideal_base_product.get_a() == ZZ(6));
    REQUIRE(quad_ideal_base_product.get_b() == ZZ(5));
    REQUIRE(quad_ideal_base_product.get_c() == ZZ(-7));

}
#endif

