#ifndef COMMON_TEST
#define COMMON_TEST

#include "catch.hpp"
#include <ANTL/common.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("Common: Kronecker Symbol computes correctly", "[Common]") {

  REQUIRE(Kronecker(long(17), long(13)) == 1);
  REQUIRE(Kronecker(long(16), long(13)) == 1);
  REQUIRE(Kronecker(long(11), long(13)) == -1);
  REQUIRE(Kronecker(long(13), long(13)) == 0);
  REQUIRE(Kronecker(long(10), long(12)) == 0);

  REQUIRE(Kronecker(ZZ(17), ZZ(13)) == 1);
  REQUIRE(Kronecker(ZZ(16), ZZ(13)) == 1);
  REQUIRE(Kronecker(ZZ(11), ZZ(13)) == -1);
  REQUIRE(Kronecker(ZZ(13), ZZ(13)) == 0);
  REQUIRE(Kronecker(ZZ(10), ZZ(12)) == 0);

  REQUIRE(Kronecker(ZZ(17), long(13)) == 1);
  REQUIRE(Kronecker(ZZ(16), long(13)) == 1);
  REQUIRE(Kronecker(ZZ(11), long(13)) == -1);
  REQUIRE(Kronecker(ZZ(13), long(13)) == 0);
  REQUIRE(Kronecker(ZZ(10), long(12)) == 0);
}

#endif
