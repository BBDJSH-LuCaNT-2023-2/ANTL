#ifndef QUADRATICORDER_ZZ_TEST
#define QUADRATICORDER_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Regulator/RegulatorLenstraData.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("RegulatorLenstra<ZZ>: Does it work?", "[RegulatorLenstra]") {

  QuadraticOrder<ZZ> quad_order{ZZ(13)};

  L_function<ZZ> l_function;
  l_function.init(ZZ(12157), 2);

  RegulatorLenstraData<ZZ> regulator_lenstra_data{&quad_order, &l_function};

  regulator_lenstra_data.regulator_lenstra();


}

#endif
