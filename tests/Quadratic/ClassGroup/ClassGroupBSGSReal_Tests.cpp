#ifndef CLASSGROUP_BSGS_TEST
#define CLASSGROUP_BSGS_TEST

#include "../../catch.hpp"
#include "../../Quadratic/TestData/TestData.hpp"

#include <ANTL/Quadratic/ClassGroup/ClassGroupBSGSReal.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>

using namespace NTL;
using namespace ANTL;

bool DBG_CGBSGSR_TEST = true;

TEST_CASE("ClassGroupReal<ZZ>: Does it work?", "[ClassGroupReal]") {

  extern const std::array<long, 1000> discriminants;
  extern const std::array<double, 1000> correct_regulators;
  extern const std::vector<std::vector<long>> correct_class_groups;

  //   int discriminants[10] = {55661, 63361, 72673, 86341, 38593,
  //                            54269, 41513, 74021, 45677, 17909};
  //
  //   double correct_regulators[10] = {25.4538649123, 17.5077374033,
  //   357.702597578,
  //                                    148.556426268,
  //                                    105.442818369, 49.458153743,
  //                                    130.11243841,  62.8055127479, 54.4615644632,
  //                                    44.1429273524};
  //
  //   vector<vector<ZZ>> correct_class_groups = {
  //       {ZZ(3)}, {ZZ(33)}, {ZZ(1)}, {ZZ(1)}, {ZZ(3)},
  //       {ZZ(1)}, {ZZ(1)},  {ZZ(1)}, {ZZ(1)}, {ZZ(1)}};

  vector<vector<long>> computed_class_groups;

  std::array<bool, 1000> computed_correctly;
  int correct_count = 0;

  for (int i = 0; i < 1000; i++) {
    // Setting up neccessary objects
    QuadraticOrder<ZZ> quad_order{ZZ(discriminants.at(i))};
    QuadraticNumber<ZZ> quad_number1{quad_order};
    QuadraticNumber<ZZ> quad_number2{quad_order};

    MultiplyComp<ZZ> mul_comp_object{};
    mul_comp_object.set_RelativeGenerator(quad_number1);
    quad_order.set_mul_comp(mul_comp_object);

    ReducePlainReal<ZZ> red_plain_real_object{};
    red_plain_real_object.set_RelativeGenerator(quad_number2);
    quad_order.set_red_best(red_plain_real_object);

    L_function<ZZ> l_function;
    l_function.init(ZZ(discriminants.at(i)), 2);

    RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order, &l_function};

    // Computing h*
    double regulator = correct_regulators.at(i);
    ZZ h_star = regulator_lenstra_data.lower_bound_hR() / regulator;

    // Setting up the ClassGroupBSGSReal object
    ClassGroupBSGSReal<ZZ> class_group_bsgs_real1{&quad_order};
    class_group_bsgs_real1.set_regulator(regulator);

    // Computing the class group
    class_group_bsgs_real1.cg_bsgs_real(h_star);

    // Adding computed class group to reslults vector
    vector<ZZ> class_group_ZZ = class_group_bsgs_real1.get_class_group();
    vector<long> class_group_long = {};
    for(auto num : class_group_ZZ) {
      class_group_long.push_back(to<long>(num));
    }
    std::sort(class_group_long.begin(), class_group_long.end());
    computed_class_groups.push_back(class_group_long);


    // Checking for corect output
    if (correct_class_groups.at(i) == computed_class_groups.at(i)) {
      computed_correctly.at(i) = true;
    } else {
      computed_correctly.at(i) = false;
    }
  }

  bool test_bool = true;

  if (DBG_CGBSGSR_TEST) {
    std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(10) << "CORRECT"
              << std::setw(12) << "COMPUTED" << std::setw(8) << "DELTA"
              << std::endl;
    std::cout << std::setw(42) << std::setfill('=') << "" << std::endl;
  }

  for (int i = 0; i < 1000; i++) {
    if (DBG_CGBSGSR_TEST) {
      std::cout << std::setfill('0') << std::setw(4) << i + 1
                << std::setfill(' ') << std::setw(6) << computed_correctly.at(i)
                << std::setw(12) << correct_class_groups.at(i)
                << std::setw(12) << computed_class_groups.at(i)
                << std::setw(8) << discriminants.at(i)
                << std::endl;
      if(computed_correctly[i]) {
        correct_count++;
      }
    }

    test_bool = test_bool && computed_correctly.at(i);
  }

  if (DBG_CGBSGSR_TEST) {
    std::cout << correct_count << "/1000 tests passed!" << std::endl;
  }


  REQUIRE(test_bool == true);
}

#endif
