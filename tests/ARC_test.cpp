#ifndef ARC_TEST
#define ARC_TEST

#include "catch.hpp"

#include <NTL/ZZ.h>

#include <fstream>
#include <iomanip>
#include <chrono>
#include <set>

using namespace NTL;

bool DBG_LENSTRA_TEST = true;

TEST_CASE("ARC_TEST: Does it work?", "[ARC]") {

  NTL::ZZ a = ZZ(123);
  std::cout << a << std::endl;

  REQUIRE(true);
}

#endif
