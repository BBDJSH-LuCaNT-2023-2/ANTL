#ifndef INDCALC_TEST
#define INDCALC_TEST

#include <string>
#include <map>
#include <iostream>
#include <NTL/ZZ.h>
#include "../catch.hpp"
#include "ANTL/IndexCalculus/IndCalc/QuadIndCalc.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Constants.hpp"

using namespace Constants;
using namespace NTL;
using namespace ANTL;

namespace UT {
  typedef int ord_t; // fake type for the QuadraticOrder
  typedef int unit_t; // fake type for number field units
  typedef int reg_t; // fake type for the regulator

// use a mock QuadraticOrder so we dont need to introduce the dependencies from ANTL/Quadratic/QuadraticOrder
  template<class T>
  class QuadraticOrder : public IOrder {
  public:
    QuadraticOrder(IOrder const &order) {}
  };
}

std::map<std::string, std::string> get_params(std::string max_num_tests_str) {
  std::map<std::string, std::string> params {{num_relations, "2"}, {size_fb, "3"}, {bound_fb, "5"},
                                             {max_num_tests, max_num_tests_str}};
  return params;
}

/*** removing this test temporarily so I can investigate a segfault in it ***/
/*** TESTS ***/
TEST_CASE("Parameters are read from the params map", "[capturing]") {
  IOrder order = IOrder();
  long expected_num_relations = 2;
  long expected_size_fb = 3;
  long expected_bound_fb = 5;
  long expected_max_num_tests = 7;
  QuadIndCalc<UT::unit_t , UT::reg_t> ind_calc1 =
    QuadIndCalc<UT::unit_t, UT::reg_t>::create(order, get_params("7"));
  REQUIRE(ind_calc1.get_relation_generator()->get_size_fb() == expected_size_fb);
  REQUIRE(ind_calc1.get_relation_generator()->get_max_num_tests() == expected_max_num_tests);
  REQUIRE(ind_calc1.get_factor_base()->get_size_fb() == expected_size_fb);
  REQUIRE(ind_calc1.get_factor_base()->get_bound() == expected_bound_fb);
}

#endif
