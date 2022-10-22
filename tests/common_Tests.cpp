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


TEST_CASE("Common: Two mpz_z are equivalent"){
  mpz_t a,b;
  mpz_init(a);
  mpz_init(b);
  mpz_set_si(a, 100);
  mpz_set_si(b, 100);
  REQUIRE(mpz_cmp(a,b)==0);
}

TEST_CASE("Common: ZZToMpz works correctly", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_t expected_a;
  mpz_set_si(expected_a, 5);
  ZZ A = ZZ(5);
  ZZToMpz(A, a);
  
  REQUIRE(mpz_cmp(a,expected_a)==0);
}

#endif
