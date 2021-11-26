#ifndef QUADRATICNUMBER_ZZ_TEST
#define QUADRATICNUMBER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

using namespace NTL;
using namespace ANTL;


TEST_CASE("QuadraticInfElement<ZZ>: baby_step correctly traverse through the principal cycle", "[QuadraticInfElement]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(193));

    QuadraticInfElement<ZZ> quad_inf_element1 = QuadraticInfElement<ZZ>(quad_order1);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -1);

    //REQUIRE(quad_inf_element1.get_distance() == 2.598698174);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    //REQUIRE(quad_inf_element1.get_distance() == 3.32835583);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
    //REQUIRE(quad_inf_element1.get_distance() == 4.82844171);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
    //REQUIRE(quad_inf_element1.get_distance() == 6.65671165);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
    //REQUIRE(quad_inf_element1.get_distance() == 6.80572746);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
    //REQUIRE(quad_inf_element1.get_distance() == 7.85709282);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
    //REQUIRE(quad_inf_element1.get_distance() == 8.15679754);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    //REQUIRE(quad_inf_element1.get_distance() == 8.71127845);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    //REQUIRE(quad_inf_element1.get_distance() == 9.16513386);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
    //REQUIRE(quad_inf_element1.get_distance() == 9.65688343);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
    //REQUIRE(quad_inf_element1.get_distance() == 10.61682945);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
    //REQUIRE(quad_inf_element1.get_distance() == 10.94102199);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
    //REQUIRE(quad_inf_element1.get_distance() == 12.84657298);

    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
    //REQUIRE(quad_inf_element1.get_distance() == 14.26937782);
    quad_inf_element1.baby_step();

    REQUIRE(quad_inf_element1.get_qib().get_a() == 1);
    REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
    REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
    //REQUIRE(quad_inf_element1.get_distance() == 15.07631652);

}

#endif

