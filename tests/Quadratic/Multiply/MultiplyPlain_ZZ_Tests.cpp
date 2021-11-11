#ifndef MULTIPLY_PLAIN_ZZ_TEST
#define MULTIPLY_PLAIN_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Multiply/MultiplyPlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("MultiplyPlain<ZZ>: Multiplying two ideals", "[MultiplyPlain]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    MultiplyPlain<ZZ> mul_plain_object = MultiplyPlain<ZZ>();
    quad_order1.set_mul_plain(mul_plain_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base_product = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(-9), ZZ(1));
    quad_ideal_base2.assign(ZZ(17), ZZ(-9), ZZ(1));


    mul(quad_ideal_base_product, quad_ideal_base1, quad_ideal_base2);

    REQUIRE(quad_ideal_base_product.get_a() == ZZ(289));
    REQUIRE(quad_ideal_base_product.get_b() == ZZ(59) );
    REQUIRE(quad_ideal_base_product.get_c() == ZZ(3)  );

}
#endif
