#ifndef REGULATOR_LENSTRA_ZZ_TEST
#define REGULATOR_LENSTRA_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Regulator/RegulatorLenstra.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>

using namespace NTL;
using namespace ANTL;

bool DBG_LENSTRA_TEST = true;

TEST_CASE("RegulatorLenstra<ZZ>: Does it work?", "[RegulatorLenstra]") {

  int discriminants[10] = {55661, 63361, 72673, 86341, 38593,
                           54269, 41513, 74021, 45677, 17909};

  double correct_regulators[10] = {25.4538649123, 17.5077374033, 357.702597578,
                                   148.556426268, 105.442818369, 49.458153743,
                                   130.11243841,  62.8055127479, 54.4615644632,
                                   44.1429273524};

  double computed_regulators[10];
  bool computed_correctly[10];

  for (int i = 0; i < 10; i++) {

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

    regulator_lenstra_data.regulator_lenstra();

    computed_regulators[i] = regulator_lenstra_data.get_regulator();

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
              << std::endl;
    std::cout << std::setw(38) << std::setfill('=') << "" << std::endl;
  }

  for (int i = 0; i < 10; i++) {
    if (DBG_LENSTRA_TEST) {
      std::cout << std::setfill('0') << std::setw(4) << i + 1
                << std::setfill(' ') << std::setw(6) << computed_correctly[i]
                << std::setw(11) << correct_regulators[i] << std::setw(10)
                << computed_regulators[i] << std::setw(7) << discriminants[i]
                << std::endl;
    }

    test_bool = test_bool && computed_correctly[i];
  }

  REQUIRE(test_bool == true);
}

#endif
