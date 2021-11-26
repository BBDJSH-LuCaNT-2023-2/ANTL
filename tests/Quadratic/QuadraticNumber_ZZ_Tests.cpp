#ifndef QUADRATICNUMBER_ZZ_TEST
#define QUADRATICNUMBER_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/Quadratic/QuadraticNumber.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("QuadraticNumber<ZZ>: Numbers are equal iff delta, a, b, d are equal", "[QuadraticNumber]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    QuadraticNumber<ZZ> quad_number1 = QuadraticNumber<ZZ>(quad_order1);
    QuadraticNumber<ZZ> quad_number2 = QuadraticNumber<ZZ>(quad_order1);

    REQUIRE(quad_number1 == quad_number2);

}

#endif
