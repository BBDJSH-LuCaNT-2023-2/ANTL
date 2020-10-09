#ifndef INDCALC_TEST
#define INDCALC_TEST

#include "ANTL/IndexCalculus/IndCalc/QuadIndCalc.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Constants.hpp"
#include <NTL/ZZ.h>
#include <string>
#include <map>
#include <iostream>

#include "../UnitTests.hpp"
#include <iostream>
#include <exception>

using namespace Constants;
using namespace NTL;

TEST_CASE("Constructors and ZZ,RR type tests", "[capturing]"){
  IOrder order = IOrder();
  long expected_num_relations = 2;
  long expected_size_fb = 3;
  long expected_bound_fb = 5;
  std::map<std::string, std::string> params {{num_relations, "2"}, {size_fb, "3"}, {bound_fb, "5"}};

  QuadIndCalc <ZZ, ZZ> ind_calc = QuadIndCalc<ZZ, ZZ>(order, params);
  REQUIRE(ind_calc.factor_base.get_size_fb() == expected_size_fb);
  REQUIRE(ind_calc.factor_base.get_bound() == expected_bound_fb);
  REQUIRE(ind_calc.relation_generator.get_size_fb() == expected_size_fb);
}

#endif