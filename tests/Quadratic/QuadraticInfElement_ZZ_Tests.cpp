#ifndef QUADRATICNUMBER_ZZ_TEST
#define QUADRATICNUMBER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticInfElement<ZZ, double>: baby_step correctly traverses the "
          "principal cycle",
          "[QuadraticInfElement]") {

  QuadraticOrder<ZZ> quad_order1{ZZ(193)};
  QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -1);

  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 2.598698174) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 3.32835583) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 4.82844171) < 0.00000001);
  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 6.65671165) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 6.80572746) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 7.85709282) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 8.15679754) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 8.71127845) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 9.16513386) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 9.65688343) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 10.61682945) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 10.94102199) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 12.84657298) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 14.26937782) < 0.00000001);

  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 1);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 15.07631652) < 0.00000001);
}

TEST_CASE("QuadraticInfElement<ZZ, double>: giant_step correctly traverses the "
          "principal cycle",
          "[QuadraticInfElement]") {

  QuadraticOrder<ZZ> quad_order1{ZZ(193)};
  QuadraticNumber<ZZ> quad_number1{quad_order1};
  QuadraticNumber<ZZ> quad_number2{quad_order1};

  MultiplyNucomp<ZZ> mul_nucomp_object{};
  mul_nucomp_object.set_RelativeGenerator(quad_number1);
  quad_order1.set_mul_nucomp(mul_nucomp_object);

  ReducePlainReal<ZZ> red_plain_real_object{};
  red_plain_real_object.set_RelativeGenerator(quad_number2);
  quad_order1.set_red_best(red_plain_real_object);

  QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};
  QuadraticInfElement<ZZ, double> quad_inf_element2{quad_order1};

  // starting point for giant steps
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();

  // ideal to be used for taking giant steps
  quad_inf_element2.baby_step();
  quad_inf_element2.baby_step();
  quad_inf_element2.baby_step();

  quad_inf_element1.giant_step(quad_inf_element2);

  REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 9.65688343) < 0.00000001);

  quad_inf_element1.giant_step(quad_inf_element2);

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 14.26937782) < 0.00000001);

  quad_inf_element1.giant_step(quad_inf_element2);

  REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 18.40467235) < 0.00000001);
}

TEST_CASE(
    "QuadraticInfElement<ZZ, double>: inverse_rho() correctly traverses the "
    "principal cycle backwards",
    "[QuadraticInfElement]") {

  QuadraticOrder<ZZ> quad_order1{ZZ(193)};
  QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};

  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();
  quad_inf_element1.baby_step();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 1);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 15.07631652) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 14.26937782) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 12.84657298) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 10.94102199) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 10.61682945) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 9.65688343) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 9.16513386) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 8.71127845) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 5);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -7);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 8.15679754) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 7.85709282) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 4);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -9);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 6.80572746) < 0.00000001);
  ;

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -2);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 6.65671165) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 2);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -3);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 4.82844171) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 3.32835583) < 0.00000001);

  quad_inf_element1.inverse_rho();

  REQUIRE(quad_inf_element1.get_qib().get_a() == 6);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 13);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -1);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 2.598698174) < 0.00000001);
}

TEST_CASE("QuadraticInfElement<ZZ, double>: adjust correctly adjusts the "
          "principal cycle",
          "[QuadraticInfElement]") {

  QuadraticOrder<ZZ> quad_order1{ZZ(193)};
  QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};

  quad_inf_element1.adjust(ZZ(4));

  REQUIRE(quad_inf_element1.get_qib().get_a() == 3);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 11);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -6);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 3.32835583) < 0.00000001);

  quad_inf_element1.adjust(ZZ(8));

  REQUIRE(quad_inf_element1.get_qib().get_a() == 7);
  REQUIRE(quad_inf_element1.get_qib().get_b() == 9);
  REQUIRE(quad_inf_element1.get_qib().get_c() == -4);
  REQUIRE(ANTL::abs(quad_inf_element1.get_distance() - 7.85709282) < 0.00000001);
}

#endif
