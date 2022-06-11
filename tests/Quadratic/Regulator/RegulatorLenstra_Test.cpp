#ifndef REGULATOR_LENSTRA_ZZ_TEST
#define REGULATOR_LENSTRA_ZZ_TEST

#include "../../Quadratic/TestData/TestData.hpp"
#include "../../catch.hpp"

#include <ANTL/Quadratic/Regulator/RegulatorLenstra.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>

using namespace NTL;
using namespace ANTL;

bool DBG_LENSTRA_TEST = true;

TEST_CASE("RegulatorLenstra<ZZ>: Does it work?", "[RegulatorLenstra]") {

  extern const std::array<long, 1100> discriminants;
  extern const std::array<double, 1100> correct_regulators;

  QuadraticOrder<ZZ> quad_order1{ZZ(30061)};
  QuadraticInfElement<ZZ, double> quad_inf_element1{quad_order1};

  int correct_count = 0;
  int test_start = 0;
  int test_bound = correct_regulators.size();

  double computed_regulators[test_bound];
  bool computed_correctly[test_bound];
  std::vector<std::string> case_types{size_t(test_bound), ""};

  for (int i = test_start; i < test_bound; i++) {
    std::cout << "Running Lenstra test " << i << " expected regulator is " << correct_regulators[i] << std::endl;
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

    QuadraticInfElement<ZZ, double> qif_1{quad_order};

    L_function<ZZ> l_function;
    l_function.init(ZZ(discriminants[i]), 2);

    RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order,
                                                            &l_function};

    ZZ bound = ZZ(0);
    regulator_lenstra_data.regulator_lenstra();
//     regulator_lenstra_data.regulator_bsgs(bound);

    computed_regulators[i] = regulator_lenstra_data.get_regulator();
    case_types.at(i) = regulator_lenstra_data.get_case_type();

    if (std::abs(computed_regulators[i] - correct_regulators[i]) < 0.000001) {
      computed_correctly[i] = true;
    } else {
      computed_correctly[i] = false;
    }
  }

  bool test_bool = true;

  if (DBG_LENSTRA_TEST) {
    std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(9) << "CORRECT"
              << std::setw(10) << "COMPUTED" << std::setw(7) << "DELTA"
              << std::setw(11) << "CASE TYPE" << std::endl;
    std::cout << std::setw(49) << std::setfill('=') << "" << std::endl;
  }

//     for(auto i : test_cases) {
  for (int i = test_start; i < test_bound; i++) {
    if (DBG_LENSTRA_TEST) {
      std::cout << std::setfill('0') << std::setw(3) << i << std::setfill(' ')
                << std::setw(6) << computed_correctly[i] << std::setw(11)
                << correct_regulators[i] << std::setw(10)
                << computed_regulators[i] << std::setw(7) << discriminants[i]
                << std::setw(35) << case_types.at(i) << std::endl;

      if (computed_correctly[i]) {
        correct_count++;
      }
    }

    test_bool = test_bool && computed_correctly[i];
  }

  if (DBG_LENSTRA_TEST) {
    std::cout << correct_count << "/" << test_bound - test_start
              << " tests passed!" << std::endl;
  }

  REQUIRE(test_bool == true);
}

#endif
