#ifndef MULTIPLY_PLAIN_LONG_TEST
#define MULTIPLY_PLAIN_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Multiply/MultiplyPlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("MultiplyPlain<long>: Multiplying two ideals", "[MultiplyPlain]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));

    MultiplyPlain<long> mul_plain_object = MultiplyPlain<long>();
    quad_order1.set_mul_plain(mul_plain_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base2 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base_product = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, -9, 1);
    quad_ideal_base2.assign(17, -9, 1);


    quad_order1.get_mul_plain()->multiply(quad_ideal_base_product, quad_ideal_base1, quad_ideal_base2);

    REQUIRE(quad_ideal_base_product.get_a() == long(289));
    REQUIRE(quad_ideal_base_product.get_b() == long(59) );
    REQUIRE(quad_ideal_base_product.get_c() == long(3)  );

}
#endif
