#ifndef CLASSGROUP_BS_TEST
#define CLASSGROUP_BS_TEST

#include "../../catch.hpp"
#include "../../Quadratic/TestData/TestData.hpp"

#include <ANTL/Quadratic/ClassGroup/ClassGroupBSReal.hpp>

#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_long.hpp>

#include <fstream>
#include <sstream>

using namespace NTL;
using namespace ANTL;

bool DBG_CGBSGSR_TEST = true;

TEST_CASE("ClassGroupReal<ZZ>: Does it work?", "[ClassGroupReal][ZZ]") {

  extern const std::vector<long> discriminants;
  extern const std::vector<double> correct_regulators;
  extern const std::vector<std::vector<long>> correct_class_groups;

  int test_method = 1;

  if(test_method == 0) {
    int correct_count = 0;
    int test_start = 0;
    int test_bound = correct_class_groups.size();
    //   int test_bound = 12;

    vector<vector<long>> computed_class_groups;

    std::vector<bool> computed_correctly;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = test_start; i < test_bound; i++) {
      std::cout << "computing test " << i << " discriminant is " << discriminants.at(i) << std::endl;
      // Setting up neccessary objects
      QuadraticOrder<ZZ> quad_order{ZZ(discriminants.at(i))};
      QuadraticNumber<ZZ> quad_number1{quad_order};
      QuadraticNumber<ZZ> quad_number2{quad_order};
      QuadraticNumber<ZZ> quad_number3{quad_order};

//       MultiplyComp<ZZ> mul_comp_object{};
//       mul_comp_object.set_RelativeGenerator(quad_number1);
//       quad_order.set_mul_comp(mul_comp_object);
//
//       ReducePlainReal<ZZ> red_plain_real_object{};
//       red_plain_real_object.set_RelativeGenerator(quad_number2);
//       quad_order.set_red_best(red_plain_real_object);

      MultiplyNucompOpt<ZZ> mul_nucomp_opt_object{};
      mul_nucomp_opt_object.set_RelativeGenerator(quad_number1);
      quad_order.set_mul_nucomp_opt(mul_nucomp_opt_object);

      ReducePlainRealOpt<ZZ> red_plain_real_opt_object{};
      red_plain_real_opt_object.set_RelativeGenerator(quad_number2);
      quad_order.set_red_best(red_plain_real_opt_object);

      SquareNuduplOpt<ZZ> sqr_nudupl_opt_object{};
      sqr_nudupl_opt_object.set_RelativeGenerator(quad_number3);
      quad_order.set_sqr_best(sqr_nudupl_opt_object);

      L_function<ZZ> l_function;
      l_function.init(ZZ(discriminants.at(i)), 2);

      RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order, &l_function};

      // Computing h*
      double regulator = correct_regulators.at(i);
      RR h_star_close = to_RR(regulator_lenstra_data.lower_bound_hR()) / to_RR(regulator);
      ZZ h_star = CeilToZZ(h_star_close);

      // Setting up the ClassGroupBSGSReal object
      ClassGroupBSReal<ZZ> class_group_bsgs_real1{&quad_order};
      class_group_bsgs_real1.set_regulator(regulator);

      // Computing the class group
      class_group_bsgs_real1.cg_bs_real(h_star);

      // Adding computed class group to reslults vector
      vector<ZZ> class_group_ZZ = class_group_bsgs_real1.get_class_group();
      vector<long> class_group_long = {};

      for(auto num : class_group_ZZ) {
        class_group_long.push_back(to<long>(num));
      }

      std::sort(class_group_long.begin(), class_group_long.end());
      std:cout << class_group_long << std::endl;
      computed_class_groups.push_back(class_group_long);

      // Checking for corect output
      if (correct_class_groups.at(i) == computed_class_groups.at(i-test_start)) {
        computed_correctly.push_back(true);
      } else {
        computed_correctly.push_back(false);
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    bool test_bool = true;

    if (DBG_CGBSGSR_TEST) {
      std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(10) << "CORRECT"
                << std::setw(21) << "COMPUTED" << std::setw(21) << "DELTA"
                << std::endl;
      std::cout << std::setw(65) << std::setfill('=') << "" << std::endl;
    }

    for (int i = test_start; i < test_bound; i++) {
      if (DBG_CGBSGSR_TEST) {
        std::cout << std::setfill('0') << std::setw(4) << i + 1
                  << std::setfill(' ') << std::setw(6) << computed_correctly.at(i)
                  << std::setw(21) << correct_class_groups.at(i)
                  << std::setw(21) << computed_class_groups.at(i)
                  << std::setw(14) << discriminants.at(i)
                  << std::endl;
        if(computed_correctly[i]) {
          correct_count++;
        }
      }

      test_bool = test_bool && computed_correctly.at(i);
    }

    if (DBG_CGBSGSR_TEST) {
      std::cout << correct_count << "/" << test_bound - test_start << " tests passed!" << std::endl;
      std::cout << "Time taken: " << duration.count() <<  " ms" << endl;
      for(int i = test_start; i < test_bound; i++) {
        if(!computed_correctly.at(i)){
          std::cout << "case " << i << " was wrong!" << std::endl;
        }
      }
    }


    REQUIRE(test_bool == true);
  }

  if(test_method == 1) {
    std::ifstream test_data;
    test_data.open("data.txt", std::ifstream::in);

    long discriminant;
    double regulator;
    std::string correct_class_group;

    int correct_count = 0;
    int test_start = 0;
    int test_bound = 0;

    std::vector<std::string> correct_testdata_class_groups;
    std::vector<std::string> computed_class_groups;
    std::vector<bool> computed_correctly;

    auto start = std::chrono::high_resolution_clock::now();
    while(test_data >> discriminant >> regulator) {
      getline(test_data, correct_class_group);
      correct_class_group.erase(0,1);
      correct_testdata_class_groups.push_back(correct_class_group);
      if(test_bound % 10000 == 0) {
        std::cout << "Doing test cases " << test_bound << " - " << test_bound + 9999 << std::endl;
      }
      test_bound++;

      // Setting up neccessary objects
      QuadraticOrder<ZZ> quad_order{ZZ(discriminant)};
      QuadraticNumber<ZZ> quad_number1{quad_order};
      QuadraticNumber<ZZ> quad_number2{quad_order};
      QuadraticNumber<ZZ> quad_number3{quad_order};

//       MultiplyComp<ZZ> mul_comp_object{};
//       mul_comp_object.set_RelativeGenerator(quad_number1);
//       quad_order.set_mul_comp(mul_comp_object);
//
//       ReducePlainReal<ZZ> red_plain_real_object{};
//       red_plain_real_object.set_RelativeGenerator(quad_number2);
//       quad_order.set_red_best(red_plain_real_object);

      MultiplyNucompOpt<ZZ> mul_nucomp_opt_object{};
      mul_nucomp_opt_object.set_RelativeGenerator(quad_number1);
      quad_order.set_mul_nucomp_opt(mul_nucomp_opt_object);

      ReducePlainRealOpt<ZZ> red_plain_real_opt_object{};
      red_plain_real_opt_object.set_RelativeGenerator(quad_number2);
      quad_order.set_red_best(red_plain_real_opt_object);

      SquareNuduplOpt<ZZ> sqr_nudupl_opt_object{};
      sqr_nudupl_opt_object.set_RelativeGenerator(quad_number3);
      quad_order.set_sqr_best(sqr_nudupl_opt_object);

      L_function<ZZ> l_function;
      l_function.init(ZZ(discriminant), 2);

      RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order, &l_function};

      // Computing h*
      RR h_star_close = to_RR(regulator_lenstra_data.lower_bound_hR()) / to_RR(regulator);
      ZZ h_star = CeilToZZ(h_star_close);

      // Setting up the ClassGroupBSGSReal object
      ClassGroupBSReal<ZZ> class_group_bsgs_real1{&quad_order};
      class_group_bsgs_real1.set_regulator(regulator);

      // Computing the class group
      class_group_bsgs_real1.cg_bs_real(h_star);

      // Adding computed class group to reslults vector
      vector<ZZ> class_group_ZZ = class_group_bsgs_real1.get_class_group();
      vector<long> class_group_long = {};

      for(auto num : class_group_ZZ) {
        class_group_long.push_back(to<long>(num));
      }

      std::sort(class_group_long.begin(), class_group_long.end());
      std::stringstream class_group_string;
      class_group_string << class_group_long;
      computed_class_groups.push_back(class_group_string.str());

      // Checking for corect output
      if (correct_testdata_class_groups.back() == computed_class_groups.back()) {
        computed_correctly.push_back(true);
      } else {
        std::cout << "ERROR AT CASE " << test_bound - 1 << std::endl;
        std::cout << "Computed: "<< computed_class_groups.back() << std::endl;
        std::cout << "Correct:  "<< correct_testdata_class_groups.back() << std::endl;
        computed_correctly.push_back(false);
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    bool test_bool = true;

//     if (DBG_CGBSGSR_TEST) {
//       std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(10) << "CORRECT"
//                 << std::setw(21) << "COMPUTED" << std::setw(21) << "DELTA"
//                 << std::endl;
//       std::cout << std::setw(65) << std::setfill('=') << "" << std::endl;
//     }

    for (int i = test_start; i < test_bound; i++) {
//       if (DBG_CGBSGSR_TEST) {
//         std::cout << std::setfill('0') << std::setw(4) << i + 1
//                   << std::setfill(' ') << std::setw(6) << computed_correctly.at(i)
//                   << std::setw(21) << correct_class_groups.at(i)
//                   << std::setw(21) << computed_class_groups.at(i)
//                   << std::setw(14) << discriminants.at(i)
//                   << std::endl;
        if(computed_correctly[i]) {
          correct_count++;
        }
//       }

      test_bool = test_bool && computed_correctly.at(i);
    }

    if (DBG_CGBSGSR_TEST) {
      std::cout << correct_count << "/" << test_bound - test_start << " tests passed!" << std::endl;
      std::cout << "Time taken: " << duration.count() <<  " ms" << endl;
      for(int i = test_start; i < test_bound; i++) {
        if(!computed_correctly.at(i)){
          std::cout << "case " << i << " was wrong!" << std::endl;
        }
      }
    }
    REQUIRE(test_bool == true);
  }

}

TEST_CASE("ClassGroupReal<long>: Does it work?", "[ClassGroupReal][long]") {

  extern const std::vector<long> discriminants;
  extern const std::vector<double> correct_regulators;
  extern const std::vector<std::vector<long>> correct_class_groups;

  int test_method = 1;

  if(test_method == 0) {
    int correct_count = 0;
    int test_start = 0;
    int test_bound = correct_class_groups.size();
    //   int test_bound = 12;

    vector<vector<long>> computed_class_groups;

    std::vector<bool> computed_correctly;

    auto start = std::chrono::high_resolution_clock::now();
    for (int i = test_start; i < test_bound; i++) {
      std::cout << "computing test " << i << " discriminant is " << discriminants.at(i) << std::endl;
      // Setting up neccessary objects
      QuadraticOrder<long> quad_order{long(discriminants.at(i))};
      QuadraticNumber<long> quad_number1{quad_order};
      QuadraticNumber<long> quad_number2{quad_order};
      QuadraticNumber<long> quad_number3{quad_order};

//       MultiplyComp<long> mul_comp_object{};
//       mul_comp_object.set_RelativeGenerator(quad_number1);
//       quad_order.set_mul_comp(mul_comp_object);
//
//       ReducePlainReal<long> red_plain_real_object{};
//       red_plain_real_object.set_RelativeGenerator(quad_number2);
//       quad_order.set_red_best(red_plain_real_object);

      MultiplyNucompOpt<long> mul_nucomp_opt_object{};
      mul_nucomp_opt_object.set_RelativeGenerator(quad_number1);
      quad_order.set_mul_nucomp_opt(mul_nucomp_opt_object);

      ReducePlainRealOpt<long> red_plain_real_opt_object{};
      red_plain_real_opt_object.set_RelativeGenerator(quad_number2);
      quad_order.set_red_best(red_plain_real_opt_object);

      SquareNuduplOpt<long> sqr_nudupl_opt_object{};
      sqr_nudupl_opt_object.set_RelativeGenerator(quad_number3);
      quad_order.set_sqr_best(sqr_nudupl_opt_object);

      L_function<long> l_function;
      l_function.init(long(discriminants.at(i)), 2);

      RegulatorLenstraData<long, double> regulator_lenstra_data{&quad_order, &l_function};

      // Computing h*
      double regulator = correct_regulators.at(i);
      RR h_star_close = to_RR(regulator_lenstra_data.lower_bound_hR()) / to_RR(regulator);
      ZZ h_star = CeilToZZ(h_star_close);

      // Setting up the ClassGroupBSGSReal object
      ClassGroupBSReal<long> class_group_bsgs_real1{&quad_order};
      class_group_bsgs_real1.set_regulator(regulator);

      // Computing the class group
      class_group_bsgs_real1.cg_bs_real(h_star);

      // Adding computed class group to results vector
      vector<ZZ> class_group_ZZ = class_group_bsgs_real1.get_class_group();
      vector<long> class_group_long = {};

      for(auto num : class_group_ZZ) {
        class_group_long.push_back(to<long>(num));
      }

      std::sort(class_group_long.begin(), class_group_long.end());
      std:cout << class_group_long << std::endl;
      computed_class_groups.push_back(class_group_long);

      // Checking for corect output
      if (correct_class_groups.at(i) == computed_class_groups.at(i-test_start)) {
        computed_correctly.push_back(true);
      } else {
        computed_correctly.push_back(false);
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    bool test_bool = true;

    if (DBG_CGBSGSR_TEST) {
      std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(10) << "CORRECT"
                << std::setw(21) << "COMPUTED" << std::setw(21) << "DELTA"
                << std::endl;
      std::cout << std::setw(65) << std::setfill('=') << "" << std::endl;
    }

    for (int i = test_start; i < test_bound; i++) {
      if (DBG_CGBSGSR_TEST) {
        std::cout << std::setfill('0') << std::setw(4) << i + 1
                  << std::setfill(' ') << std::setw(6) << computed_correctly.at(i)
                  << std::setw(21) << correct_class_groups.at(i)
                  << std::setw(21) << computed_class_groups.at(i)
                  << std::setw(14) << discriminants.at(i)
                  << std::endl;
        if(computed_correctly[i]) {
          correct_count++;
        }
      }

      test_bool = test_bool && computed_correctly.at(i);
    }

    if (DBG_CGBSGSR_TEST) {
      std::cout << correct_count << "/" << test_bound - test_start << " tests passed!" << std::endl;
      std::cout << "Time taken: " << duration.count() <<  " ms" << endl;
      for(int i = test_start; i < test_bound; i++) {
        if(!computed_correctly.at(i)){
          std::cout << "case " << i << " was wrong!" << std::endl;
        }
      }
    }


    REQUIRE(test_bool == true);
  }

  if(test_method == 1) {
    std::ifstream test_data;
    test_data.open("data.txt", std::ifstream::in);

    long discriminant;
    double regulator;
    std::string correct_class_group;

    int correct_count = 0;
    int test_start = 0;
    int test_bound = 0;

    std::vector<std::string> correct_testdata_class_groups;
    std::vector<std::string> computed_class_groups;
    std::vector<bool> computed_correctly;

    auto start = std::chrono::high_resolution_clock::now();
    while(test_data >> discriminant >> regulator) {
      getline(test_data, correct_class_group);
      correct_class_group.erase(0,1);
      correct_testdata_class_groups.push_back(correct_class_group);
      if(test_bound % 10000 == 0) {
        std::cout << "Doing test cases " << test_bound << " - " << test_bound + 9999 << std::endl;
      }
      test_bound++;

      // Setting up neccessary objects
      QuadraticOrder<long> quad_order{long(discriminant)};
      QuadraticNumber<long> quad_number1{quad_order};
      QuadraticNumber<long> quad_number2{quad_order};
      QuadraticNumber<long> quad_number3{quad_order};

//       MultiplyComp<long> mul_comp_object{};
//       mul_comp_object.set_RelativeGenerator(quad_number1);
//       quad_order.set_mul_comp(mul_comp_object);
//
//       ReducePlainReal<long> red_plain_real_object{};
//       red_plain_real_object.set_RelativeGenerator(quad_number2);
//       quad_order.set_red_best(red_plain_real_object);

      MultiplyNucompOpt<long> mul_nucomp_opt_object{};
      mul_nucomp_opt_object.set_RelativeGenerator(quad_number1);
      quad_order.set_mul_nucomp_opt(mul_nucomp_opt_object);

      ReducePlainRealOpt<long> red_plain_real_opt_object{};
      red_plain_real_opt_object.set_RelativeGenerator(quad_number2);
      quad_order.set_red_best(red_plain_real_opt_object);

      SquareNuduplOpt<long> sqr_nudupl_opt_object{};
      sqr_nudupl_opt_object.set_RelativeGenerator(quad_number3);
      quad_order.set_sqr_best(sqr_nudupl_opt_object);

      L_function<long> l_function;
      l_function.init(long(discriminant), 2);

      RegulatorLenstraData<long, double> regulator_lenstra_data{&quad_order, &l_function};

      // Computing h*
      RR h_star_close = to_RR(regulator_lenstra_data.lower_bound_hR()) / to_RR(regulator);
      ZZ h_star = CeilToZZ(h_star_close);

      // Setting up the ClassGroupBSGSReal object
      ClassGroupBSReal<long> class_group_bsgs_real1{&quad_order};
      class_group_bsgs_real1.set_regulator(regulator);

      // Computing the class group
      class_group_bsgs_real1.cg_bs_real(h_star);

      // Adding computed class group to results vector
      vector<ZZ> class_group_ZZ = class_group_bsgs_real1.get_class_group();
      vector<long> class_group_long = {};

      for(auto num : class_group_ZZ) {
        class_group_long.push_back(to<long>(num));
      }

      std::sort(class_group_long.begin(), class_group_long.end());
      std::stringstream class_group_string;
      class_group_string << class_group_long;
      computed_class_groups.push_back(class_group_string.str());

      // Checking for corect output
      if (correct_testdata_class_groups.back() == computed_class_groups.back()) {
        computed_correctly.push_back(true);
      } else {
        std::cout << "ERROR AT CASE " << test_bound - 1 << std::endl;
        std::cout << "Computed: "<< computed_class_groups.back() << std::endl;
        std::cout << "Correct:  "<< correct_testdata_class_groups.back() << std::endl;
        computed_correctly.push_back(false);
      }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    bool test_bool = true;

//     if (DBG_CGBSGSR_TEST) {
//       std::cout << "CASE" << std::setw(8) << "RESULT" << std::setw(10) << "CORRECT"
//                 << std::setw(21) << "COMPUTED" << std::setw(21) << "DELTA"
//                 << std::endl;
//       std::cout << std::setw(65) << std::setfill('=') << "" << std::endl;
//     }

    for (int i = test_start; i < test_bound; i++) {
//       if (DBG_CGBSGSR_TEST) {
//         std::cout << std::setfill('0') << std::setw(4) << i + 1
//                   << std::setfill(' ') << std::setw(6) << computed_correctly.at(i)
//                   << std::setw(21) << correct_class_groups.at(i)
//                   << std::setw(21) << computed_class_groups.at(i)
//                   << std::setw(14) << discriminants.at(i)
//                   << std::endl;
        if(computed_correctly[i]) {
          correct_count++;
        }
//       }

      test_bool = test_bool && computed_correctly.at(i);
    }

    if (DBG_CGBSGSR_TEST) {
      std::cout << correct_count << "/" << test_bound - test_start << " tests passed!" << std::endl;
      std::cout << "Time taken: " << duration.count() <<  " ms" << endl;
      for(int i = test_start; i < test_bound; i++) {
        if(!computed_correctly.at(i)){
          std::cout << "case " << i << " was wrong!" << std::endl;
        }
      }
    }
    REQUIRE(test_bool == true);
  }

}

#endif
