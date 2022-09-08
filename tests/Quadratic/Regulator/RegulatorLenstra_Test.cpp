#ifndef REGULATOR_LENSTRA_ZZ_TEST
#define REGULATOR_LENSTRA_ZZ_TEST

#include "../../Quadratic/TestData/TestData.hpp"
#include "../../catch.hpp"

#include <ANTL/Quadratic/Regulator/RegulatorLenstra.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>

#include <fstream>
#include <chrono>

using namespace NTL;
using namespace ANTL;

bool DBG_LENSTRA_TEST = true;

TEST_CASE("RegulatorLenstra<ZZ>: Does it work?", "[RegulatorLenstra]") {

  extern const std::vector<long> discriminants;
  extern const std::vector<double> correct_regulators;

  QuadraticOrder<ZZ> quad_order1{ZZ(30061)};
  QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};

  int test_method = 1;

  if(test_method == 0) {
    int correct_count = 0;
    int test_start = 0;
    int test_bound = correct_regulators.size();

    double computed_regulators[test_bound];
    bool computed_correctly[test_bound];
    ZZ computed_hstars[test_bound];

    std::vector<std::string> case_types{size_t(test_bound), ""};

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = test_start; i < test_bound; i++) {
      std::cout << "Running Lenstra test " << i << " expected regulator is "
                << correct_regulators[i] << std::endl;
      std::cout << "Discriminat is " << discriminants[i] << std::endl;
      QuadraticOrder<ZZ> quad_order{ZZ(discriminants[i])};
      QuadraticNumber<ZZ> quad_number1{quad_order};
      QuadraticNumber<ZZ> quad_number2{quad_order};

      MultiplyComp<ZZ> mul_comp_object{};
      mul_comp_object.set_RelativeGenerator(quad_number1);
      quad_order.set_mul_comp(mul_comp_object);

      ReducePlainReal<ZZ> red_plain_real_object{};
      red_plain_real_object.set_RelativeGenerator(quad_number2);
      quad_order.set_red_best(red_plain_real_object);

      L_function<ZZ> l_function;
      l_function.init(ZZ(discriminants[i]), 2);

      RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order,
                                                              &l_function};

      ZZ bound = ZZ(0);
      regulator_lenstra_data.regulator_lenstra();
      //     regulator_lenstra_data.regulator_bsgs(bound);

      computed_regulators[i] = regulator_lenstra_data.get_regulator();
      case_types.at(i) = regulator_lenstra_data.get_case_type();

      computed_hstars[i] = regulator_lenstra_data.get_hstar();

      if (std::abs(computed_regulators[i] - correct_regulators[i]) < 0.000001) {
        computed_correctly[i] = true;
      } else {
        computed_correctly[i] = false;
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    bool test_bool = true;

    if (DBG_LENSTRA_TEST) {
      std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(9) << "CORRECT"
                << std::setw(10) << "COMPUTED" << std::setw(8) << "H_STAR" << std::setw(7) << "DELTA"
                << std::setw(11) << "CASE TYPE" << std::endl;
      std::cout << std::setw(57) << std::setfill('=') << "" << std::endl;
    }

    //     for(auto i : test_cases) {
    for (int i = test_start; i < test_bound; i++) {
      if (DBG_LENSTRA_TEST) {
        std::cout << std::setfill('0') << std::setw(3) << i << std::setfill(' ')
                  << std::setw(6) << computed_correctly[i] << std::setw(11)
                  << correct_regulators[i] << std::setw(10)
                  << computed_regulators[i] << std::setw(8) << computed_hstars[i] << std::setw(7) << discriminants[i]
                  << std::setw(35) << std::endl;

        if (computed_correctly[i]) {
          correct_count++;
        }
      }

      test_bool = test_bool && computed_correctly[i];
    }

    if (DBG_LENSTRA_TEST) {
      std::cout << correct_count << "/" << test_bound - test_start
                << " tests passed!" << std::endl;
      std::cout << "Time taken: " << duration.count() <<  " ms" << endl;
    }

    REQUIRE(test_bool == true);
  }

  if(test_method == 1) {

    std::ifstream test_data;
    test_data.open("data.txt", std::ifstream::in);

    long discriminant;
    double regulator;
    std::string class_group;

    int correct_count = 0;
    int test_start = 0;
    int test_bound = 0;

    std::vector<double> correct_testdata_regulators;
    std::vector<double> computed_regulators;
    std::vector<bool> computed_correctly;

    std::vector<std::string> case_types{size_t(test_bound), ""};

    auto start = std::chrono::high_resolution_clock::now();
    while(test_data >> discriminant >> regulator) {
      getline(test_data, class_group);
      correct_testdata_regulators.push_back(regulator);
      if(test_bound % 10000 == 0) {
        std::cout << "Doing test cases " << test_bound << " - " << test_bound + 9999 << std::endl;
      }
      test_bound++;

      QuadraticOrder<ZZ> quad_order{ZZ(discriminant)};
      QuadraticNumber<ZZ> quad_number1{quad_order};
      QuadraticNumber<ZZ> quad_number2{quad_order};

      MultiplyComp<ZZ> mul_comp_object{};
      mul_comp_object.set_RelativeGenerator(quad_number1);
      quad_order.set_mul_comp(mul_comp_object);

      ReducePlainReal<ZZ> red_plain_real_object{};
      red_plain_real_object.set_RelativeGenerator(quad_number2);
      quad_order.set_red_best(red_plain_real_object);

      L_function<ZZ> l_function;
      l_function.init(ZZ(discriminant), 2);

      RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order,
                                                              &l_function};

      ZZ bound = ZZ(0);
      regulator_lenstra_data.regulator_lenstra();
      //     regulator_lenstra_data.regulator_bsgs(bound);

      correct_testdata_regulators.push_back(regulator_lenstra_data.get_regulator());
      case_types.push_back(regulator_lenstra_data.get_case_type());

      if (std::abs(correct_testdata_regulators.back() - regulator) < 0.000001) {
        computed_correctly.push_back(true);
      } else {
        std::cout << "ERROR AT CASE " << test_bound - 1 << std::endl;
        computed_correctly.push_back(false);
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    bool test_bool = true;

//     if (DBG_LENSTRA_TEST) {
//       std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(9) << "CORRECT"
//                 << std::setw(10) << "COMPUTED" << std::setw(8)  << "DELTA"
//                 << std::setw(11) << "CASE TYPE" << std::endl;
//       std::cout << std::setw(57) << std::setfill('=') << "" << std::endl;
//     }

    //     for(auto i : test_cases) {
    for (int i = test_start; i < test_bound; i++) {
//       if (DBG_LENSTRA_TEST) {
//         std::cout << std::setfill('0') << std::setw(3) << i << std::setfill(' ')
//                   << std::setw(6) << computed_correctly[i] << std::setw(11)
//                   << correct_regulators[i] << std::setw(10)
//                   << computed_regulators[i] << std::setw(8) << discriminants[i]
//                   << std::setw(35) << std::endl;
//
        if (computed_correctly[i]) {
          correct_count++;
        }
//       }

      test_bool = test_bool && computed_correctly[i];
    }

    if (DBG_LENSTRA_TEST) {
      std::cout << correct_count << "/" << test_bound - test_start
                << " tests passed!" << std::endl;
      std::cout << "Time taken: " << duration.count() <<  " ms" << endl;
    }

    REQUIRE(test_bool == true);
  }

}

#endif
