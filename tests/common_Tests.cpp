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


TEST_CASE("Common: ZZToMpz Basic case", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_t expected_a;
  mpz_set_si(expected_a, 5);
  ZZ A = ZZ(5);
  ZZToMpz(A, a);
  
  REQUIRE(mpz_cmp(a,expected_a)==0);
}

TEST_CASE("Common: ZZToMpz Negative case", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_t expected_a;
  mpz_set_si(expected_a, -101);
  ZZ A = ZZ(-101);
  ZZToMpz(A, a);
  //cout << "A: " << A;
  //gmp_printf(" a: %Zd\n", a);
  //gmp_printf("expected_a: %Zd\n", expected_a);

  REQUIRE(mpz_cmp(a, expected_a) == 0);
}

TEST_CASE("Common: ZZToMpz Zero Case", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_t expected_a;
  mpz_set_si(expected_a, 0);
  ZZ A = ZZ(0);
  ZZToMpz(A, a);

  REQUIRE(mpz_cmp(a, expected_a) == 0);
}
TEST_CASE("Common: ZZToMpz <128 bit test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_t expected_a;
  mpz_init(expected_a);
  mpz_set_str(expected_a, "249208652633223703531132231181060005980",10);
  ZZ A;
  A = conv<ZZ>("249208652633223703531132231181060005980");
  ZZToMpz(A, a);

  REQUIRE(mpz_cmp(a, expected_a) == 0);
}
TEST_CASE("Common: ZZToMpz <256 bit test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_t expected_a;
  mpz_init(expected_a);
  mpz_set_str(expected_a, "58963287023509755963336981841255765012854373639575441005988935943903024514155", 10);
  ZZ A;
  A = conv<ZZ>("58963287023509755963336981841255765012854373639575441005988935943903024514155");
  ZZToMpz(A, a);

  REQUIRE(mpz_cmp(a, expected_a) == 0);
}

TEST_CASE("Common: MpzToZZ Basic Test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_set_si(a, 5);
  ZZ expected_A;
  expected_A = 5;
  ZZ A;
  MpzToZZ(a,A);

  REQUIRE(A == expected_A);

}

TEST_CASE("Common: MpzToZZ Negative Test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_set_si(a, -10);
  ZZ expected_A;
  expected_A = -10;
  ZZ A;
  MpzToZZ(a,A);

  REQUIRE(A == expected_A);
}

TEST_CASE("Common: MpzToZZ Zero Test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_set_si(a, 0);
  ZZ expected_A;
  expected_A = 0;
  ZZ A;
  MpzToZZ(a, A);

  REQUIRE(A == expected_A);
}

TEST_CASE("Common: MpzToZZ <128 bit test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_set_str(a, "263554072397940936318525601434069741669", 10);
  ZZ expected_A = conv<ZZ>("263554072397940936318525601434069741669");
  ZZ A;
  MpzToZZ(a, A);

  REQUIRE(A == expected_A);
}
TEST_CASE("Common: MpzToZZ <256 bit test", "[Common]"){
  mpz_t a;
  mpz_init(a);
  mpz_set_str(a, "96647711664069877376733932946779407802298220209015413512147548816784005205108", 10);
  ZZ expected_A = conv<ZZ>("96647711664069877376733932946779407802298220209015413512147548816784005205108");
  ZZ A;
  MpzToZZ(a, A);

  REQUIRE(A == expected_A);
}
#endif
