#ifndef QUADRATICORDER_ZZ_TEST
#define QUADRATICORDER_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Regulator/RegulatorLenstra.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("RegulatorLenstra<ZZ>: Does it work?", "[RegulatorLenstra]") {

  QuadraticOrder<ZZ> quad_order{ZZ(12157)};
  QuadraticNumber<ZZ> quad_number1{quad_order};
  QuadraticNumber<ZZ> quad_number2{quad_order};

  MultiplyNucomp<ZZ> mul_nucomp_object{};
  mul_nucomp_object.set_RelativeGenerator(quad_number1);
  quad_order.set_mul_nucomp(mul_nucomp_object);

  ReducePlainReal<ZZ> red_plain_real_object{};
  red_plain_real_object.set_RelativeGenerator(quad_number2);
  quad_order.set_red_best(red_plain_real_object);


  L_function<ZZ> l_function;
  l_function.init(ZZ(12157), 2);

  RegulatorLenstraData<ZZ, double> regulator_lenstra_data{&quad_order, &l_function};

  regulator_lenstra_data.regulator_lenstra();

  std::cout << "regulator is " << regulator_lenstra_data.get_regulator() << std::endl;
  std::cout << "it should be " << 43.7136859265 << std::endl;

  REQUIRE(ANTL::abs(regulator_lenstra_data.get_regulator() - 43.7136859265) < 0.00000001);
}

#endif
