#ifndef CLASSGROUP_BSGS_TEST
#define CLASSGROUP_BSGS_TEST

#include "../../catch.hpp"

#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>
#include <ANTL/Quadratic/ClassGroup/ClassGroupBSGSReal.hpp>

using namespace NTL;
using namespace ANTL;

bool DBG_CGBSGSR_TEST = false;

TEST_CASE("ClassGroupReal<ZZ>: Does it work?", "[RegulatorLenstra]") {

  int discriminants[10] = {55661, 63361, 72673, 86341, 38593,
                           54269, 41513, 74021, 45677, 17909};

  double correct_regulators[10] = {25.4538649123, 17.5077374033, 357.702597578,
                                   148.556426268, 105.442818369, 49.458153743,
                                   130.11243841,  62.8055127479, 54.4615644632,
                                   44.1429273524};

  double computed_regulators[10];
  bool computed_correctly[10];

  int i = 0;

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

  ClassGroupBSGSReal<ZZ> class_group_bsgs_real1{&quad_order};

  ZZ h_star = regulator_lenstra_data.lower_bound_hR() /
              regulator_lenstra_data.get_regulator();
  class_group_bsgs_real1.cg_bsgs_real(h_star);

  REQUIRE(1 == 1);
}

#endif
