#ifndef LFUNCTION_ZZ_TEST
#define LFUNCTION_ZZ_TEST

#include "../catch.hpp"
#include <ANTL/L_function/L_function.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("LFunction<ZZ>: Computing L(1, Chi(Delta))", "[LFunction]") {

    L_function<ZZ> L_function1;

    L_function1.init(ZZ(12157), 2);

    RR LFunctionApprox1 = L_function1.approximate(1,long(10000000));

    REQUIRE(abs(LFunctionApprox1 - 0.792928842) < 0.00001);
}

#endif

