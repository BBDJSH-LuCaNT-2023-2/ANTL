#ifndef QUADRATICNUMBER_ZZ_TEST
#define QUADRATICNUMBER_ZZ_TEST

#include "TestData/TestData.hpp"
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
/*
TEST_CASE("QuadraticInfElement<ZZ, double>: giant_step correctly traverses the "
          "principal cycle",
          "[QuadraticInfElement]") {

  QuadraticOrder<ZZ> quad_order1{ZZ(193)};
  QuadraticNumber<ZZ> quad_number1{quad_order1};
  QuadraticNumber<ZZ> quad_number2{quad_order1};

  MultiplyComp<ZZ> mul_comp_object{};
  mul_comp_object.set_RelativeGenerator(quad_number1);
  quad_order1.set_mul_comp(mul_comp_object);

  ReducePlainRealOpt<ZZ> red_plain_real_object{};
//   ReducePlainReal<ZZ> red_plain_real_object{};
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
*/

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
/*
TEST_CASE("QuadraticInfElement<ZZ, double>: Testing Nudupl_Opt and ReducePlainReal_Opt",
          "[QuadraticInfElement]") {

  extern const std::vector<long> discriminants;

  //
  std::vector<bool> case_results;

  for(int i = 0; i < discriminants.size(); i++) {

    QuadraticOrder<ZZ> quad_order1{ZZ(discriminants.at(i))};
    QuadraticNumber<ZZ> quad_number{quad_order1};

    std::cout << "Case: D = " << discriminants.at(i) << std::endl;
    SquareNuduplOpt<ZZ> sqr_nudupl_opt_object{};
    sqr_nudupl_opt_object.set_RelativeGenerator(quad_number);
    quad_order1.set_sqr_best(sqr_nudupl_opt_object);

    ReducePlainRealOpt<ZZ> red_plain_real_object{};
    red_plain_real_object.set_RelativeGenerator(quad_number);
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticInfElement<ZZ, double> quad_inf_element{quad_order1};


    // compute infrastructure (twice) and regulator
    std::cout << "baby step list:" << std::endl;
    std::vector<QuadraticInfElement<ZZ, double>> infrastructure;

    do {
      std::cout << quad_inf_element.get_qib() << ", " << quad_inf_element.get_distance() << std::endl;
      infrastructure.push_back(quad_inf_element);
      quad_inf_element.baby_step();
    } while(!quad_inf_element.is_one());

    std::cout << quad_inf_element.get_qib() << ", " << quad_inf_element.get_distance() << std::endl;
    infrastructure.push_back(quad_inf_element);
    double regulator = quad_inf_element.get_distance();

    while(quad_inf_element.get_distance() <= 2*regulator) {
      std::cout << quad_inf_element.get_qib() << ", " << quad_inf_element.get_distance() << std::endl;
      infrastructure.push_back(quad_inf_element);
      quad_inf_element.baby_step();
    }

    // repeatedly square a_2 (that is, rho(identity)) until past regulator
    // and ensure operational integrity
    quad_inf_element.assign_one();
    quad_inf_element.baby_step();
    quad_inf_element.giant_step(quad_inf_element);

    bool correct_behaviour = true;


    while(quad_inf_element.get_distance() <= regulator + 0.00000001) {
      // Check if quad_inf_element is in the infrastructure
      bool is_in = false;
      for (auto element : infrastructure) {
        if(element == quad_inf_element && ANTL::abs(element.get_distance() - quad_inf_element.get_distance()) < 0.01) {
          is_in = true;
          break;
        }
      }
      if(is_in == false) {
        std::cout << "Case " << i << " failed! Discriminant = " << discriminants.at(i) << std::endl;
        std::cout << quad_inf_element.get_qib() << ", " << quad_inf_element.get_distance() << std::endl;
        correct_behaviour = false;
        break;
      }

      // Square quad_inf_element
      quad_inf_element.giant_step(quad_inf_element);
    }

    // if correct_behaviour == false, store result and skip a_3 to continue to next case

    if(correct_behaviour == true) {
      // repeatedly square a_3 identity until past regulator
      // and ensure operational integrity
      quad_inf_element.assign_one();
      quad_inf_element.baby_step();
      quad_inf_element.baby_step();
      quad_inf_element.giant_step(quad_inf_element);


      while(quad_inf_element.get_distance() <= regulator + 0.00000001) {
        // Check if quad_inf_element is in the infrastructure
        bool is_in = false;
        for (auto element : infrastructure) {
          if(element == quad_inf_element && ANTL::abs(element.get_distance() - quad_inf_element.get_distance()) < 0.01) {
            is_in = true;
            break;
          }
        }
        if(is_in == false) {
          std::cout << "Case " << i << " failed! Discriminant = " << discriminants.at(i) << std::endl;
          std::cout << quad_inf_element.get_qib() << ", " << quad_inf_element.get_distance() << std::endl;
          correct_behaviour = false;
          break;
        }
          // Square quad_inf_element
          quad_inf_element.giant_step(quad_inf_element);
        }
    }

    // Store case result
    case_results.push_back(correct_behaviour);
  }

  bool final_result = true;
  long failed_cases = 0;

  for (auto case_result : case_results) {
    if (case_result == false) failed_cases++;
    final_result &= case_result;
  }

  std::cout << "NUDUPL_OPT cases tested: " << case_results.size() - 1 << std::endl;
  std::cout << failed_cases << "/" << discriminants.size() - 1 << " failed!" << std::endl;
  REQUIRE(final_result == true);
}*/

TEST_CASE("QuadraticInfElement<ZZ, double>: Testing Nucomp_Opt and ReducePlainReal_Opt",
          "[QuadraticInfElement]") {

  extern const std::vector<long> discriminants;

  //
  std::vector<bool> case_results;

  for(int i = 0; i < discriminants.size(); i++) {

    QuadraticOrder<ZZ> quad_order1{ZZ(discriminants.at(i))};
    QuadraticNumber<ZZ> quad_number1{quad_order1};
    QuadraticNumber<ZZ> quad_number2{quad_order1};
    QuadraticNumber<ZZ> quad_number3{quad_order1};

    std::cout << "Case: D = " << discriminants.at(i) << std::endl;

    MultiplyNucompOpt<ZZ> mul_nucomp_opt_object{};
    mul_nucomp_opt_object.set_RelativeGenerator(quad_number1);
    quad_order1.set_mul_nucomp_opt(mul_nucomp_opt_object);

    SquareNuduplOpt<ZZ> sqr_nudupl_opt_object{};
    sqr_nudupl_opt_object.set_RelativeGenerator(quad_number2);
    quad_order1.set_sqr_best(sqr_nudupl_opt_object);

    ReducePlainRealOpt<ZZ> red_plain_real_object{};
    red_plain_real_object.set_RelativeGenerator(quad_number3);
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};
    QuadraticInfElement<ZZ, double> quad_inf_element2{quad_order1};


    // compute infrastructure (twice) and regulator
    std::cout << "NUCOMP_OPT: baby step list:" << std::endl;
    std::vector<QuadraticInfElement<ZZ, double>> infrastructure;

    do {
      std::cout << quad_inf_element1.get_qib() << ", " << quad_inf_element1.get_distance() << std::endl;
      infrastructure.push_back(quad_inf_element1);
      quad_inf_element1.baby_step();
    } while(!quad_inf_element1.is_one());

    std::cout << quad_inf_element1.get_qib() << ", " << quad_inf_element1.get_distance() << std::endl;
    infrastructure.push_back(quad_inf_element1);
    double regulator = quad_inf_element1.get_distance();

    while(quad_inf_element1.get_distance() <= 2*regulator) {
      std::cout << quad_inf_element1.get_qib() << ", " << quad_inf_element1.get_distance() << std::endl;
      infrastructure.push_back(quad_inf_element1);
      quad_inf_element1.baby_step();
    }

    // repeatedly multiply by a_2 (that is, rho(identity)) until past regulator
    // and ensure operational integrity
    quad_inf_element1.assign_one();
    quad_inf_element1.baby_step();

    quad_inf_element2.assign_one();
    quad_inf_element2.giant_step(quad_inf_element1);

    std::cout << "the newly computed inf element is " << quad_inf_element2.get_qib() << ", " << quad_inf_element2.get_distance() << std::endl;

    bool correct_behaviour = true;


    while(quad_inf_element2.get_distance() <= regulator + 0.00000001) {
      // Check if quad_inf_element is in the infrastructure
      bool is_in = false;
      for (auto element : infrastructure) {
        if(element == quad_inf_element2 && ANTL::abs(element.get_distance() - quad_inf_element2.get_distance()) < 0.01) {
          is_in = true;
          break;
        }
      }
      if(is_in == false) {
        std::cout << "Case " << i << " failed! Discriminant = " << discriminants.at(i) << std::endl;
        std::cout << quad_inf_element2.get_qib() << ", " << quad_inf_element2.get_distance() << std::endl;
        correct_behaviour = false;
        break;
      }

      // Square quad_inf_element
      quad_inf_element2.giant_step(quad_inf_element1);
    }

    // if correct_behaviour == false, store result and skip a_3 to continue to next case

    if(correct_behaviour == true) {
      // repeatedly multiply by a_2 (that is, rho(identity)) until past regulator
      // and ensure operational integrity
      quad_inf_element1.assign_one();
      quad_inf_element1.baby_step();
      quad_inf_element1.baby_step();

      quad_inf_element2.assign_one();
      quad_inf_element2.giant_step(quad_inf_element1);

      bool correct_behaviour = true;

      while(quad_inf_element2.get_distance() <= regulator + 0.00000001) {
        // Check if quad_inf_element is in the infrastructure
        bool is_in = false;
        for (auto element : infrastructure) {
          if(element == quad_inf_element2 && ANTL::abs(element.get_distance() - quad_inf_element2.get_distance()) < 0.01) {
            is_in = true;
            break;
          }
        }
        if(is_in == false) {
          std::cout << "Case " << i << " failed! Discriminant = " << discriminants.at(i) << std::endl;
          std::cout << quad_inf_element2.get_qib() << ", " << quad_inf_element2.get_distance() << std::endl;
          correct_behaviour = false;
          break;
        }

        // Square quad_inf_element
        quad_inf_element2.giant_step(quad_inf_element1);
      }
    }

    // Store case result
    case_results.push_back(correct_behaviour);
  }

  bool final_result = true;
  long failed_cases = 0;

  for (auto case_result : case_results) {
    if (case_result == false) failed_cases++;
    final_result &= case_result;
  }

  std::cout << "NUCOMP_OPT cases tested: " << case_results.size() - 1 << std::endl;
  std::cout << failed_cases << "/" << discriminants.size() - 1 << " failed!" << std::endl;
  REQUIRE(final_result == true);
}

#endif
