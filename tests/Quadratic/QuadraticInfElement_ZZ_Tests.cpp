#ifndef QUADRATICNUMBER_ZZ_TEST
#define QUADRATICNUMBER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

using namespace NTL;
using namespace ANTL;


TEST_CASE("QuadraticInfElement<ZZ>: baby_step correctly traverses the principal cycle", "[QuadraticInfElement]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(193));
    QuadraticInfElement<ZZ> quad_inf_element1 = QuadraticInfElement<ZZ>(quad_order1);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -1);

    REQUIRE(abs(quad_inf_element1.get_distance() - 2.598698174) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    REQUIRE(abs(quad_inf_element1.get_distance() - 3.32835583) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
    REQUIRE(abs(quad_inf_element1.get_distance() - 4.82844171) < 0.00000001);
    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
    REQUIRE(abs(quad_inf_element1.get_distance() - 6.65671165) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
    REQUIRE(abs(quad_inf_element1.get_distance() - 6.80572746) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
    REQUIRE(abs(quad_inf_element1.get_distance() - 7.85709282) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
    REQUIRE(abs(quad_inf_element1.get_distance() - 8.15679754) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    REQUIRE(abs(quad_inf_element1.get_distance() - 8.71127845) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    REQUIRE(abs(quad_inf_element1.get_distance() - 9.16513386) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
    REQUIRE(abs(quad_inf_element1.get_distance() - 9.65688343) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
    REQUIRE(abs(quad_inf_element1.get_distance() - 10.61682945) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
    REQUIRE(abs(quad_inf_element1.get_distance() - 10.94102199) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
    REQUIRE(abs(quad_inf_element1.get_distance() - 12.84657298) < 0.00000001);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
    REQUIRE(abs(quad_inf_element1.get_distance() - 14.26937782) < 0.00000001);
    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 1);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    REQUIRE(abs(quad_inf_element1.get_distance() - 15.07631652) < 0.00000001);
}

TEST_CASE("QuadraticInfElement<ZZ>: giant_step correctly traverses the principal cycle", "[QuadraticInfElement]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(193));
    QuadraticNumber<ZZ> quad_number1 = QuadraticNumber<ZZ>(quad_order1);
    QuadraticNumber<ZZ> quad_number2 = QuadraticNumber<ZZ>(quad_order1);

    MultiplyNucomp<ZZ> mul_nucomp_object = MultiplyNucomp<ZZ>();
    mul_nucomp_object.set_RelativeGenerator(quad_number1);
    quad_order1.set_mul_nucomp(mul_nucomp_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    red_plain_real_object.set_RelativeGenerator(quad_number2);
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticInfElement<ZZ> quad_inf_element1 = QuadraticInfElement<ZZ>(quad_order1);
    QuadraticInfElement<ZZ> quad_inf_element2 = QuadraticInfElement<ZZ>(quad_order1);

    //starting point for giant steps
    quad_inf_element1.baby_step();
    quad_inf_element1.baby_step();
    quad_inf_element1.baby_step();

    //ideal to be used for taking giant steps
    quad_inf_element2.baby_step();
    quad_inf_element2.baby_step();
    quad_inf_element2.baby_step();

    quad_inf_element1.giant_step(quad_inf_element2);

    REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
    REQUIRE(abs(quad_inf_element1.get_distance() - 9.65688343) < 0.00000001);

    quad_inf_element1.giant_step(quad_inf_element2);

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
    REQUIRE(abs(quad_inf_element1.get_distance() - 14.26937782) < 0.00000001);

    quad_inf_element1.giant_step(quad_inf_element2);

    REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    REQUIRE(abs(quad_inf_element1.get_distance() - 18.40467235) < 0.00000001);
}
#endif

