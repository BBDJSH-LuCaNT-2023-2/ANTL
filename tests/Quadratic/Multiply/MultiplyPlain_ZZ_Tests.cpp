#ifndef MULTIPLY_PLAIN_ZZ_TEST
#define MULTIPLY_PLAIN_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Multiply/MultiplyPlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("MultiplyPlain<ZZ>: Multiplying two ideals", "[MultiplyPlain]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(12157));

    MultiplyPlain<ZZ> mul_plain_object = MultiplyPlain<ZZ>();
    quad_order1.set_mul_plain(mul_plain_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base2 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> qib_product = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(57), ZZ(23), ZZ(-51));
    quad_ideal_base2.assign(ZZ(51), ZZ(11), ZZ(-59));


    mul(qib_product, quad_ideal_base1, quad_ideal_base2);

    std::cout << "a, b, c are " << qib_product.get_a() << " " <<  qib_product.get_b() << " " <<  qib_product.get_c() << std::endl;
//     REQUIRE(quad_ideal_base_product.get_a() == ZZ(289));
//     REQUIRE(quad_ideal_base_product.get_b() == ZZ(59) );
//     REQUIRE(quad_ideal_base_product.get_c() == ZZ(3)  );

}
#endif
