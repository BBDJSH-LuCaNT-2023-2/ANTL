#ifndef CUBE_NUCUBE_ZZ_TEST
#define CUBE_NUCUBE_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Cube/CubeNucube.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("CubeNucube<ZZ>: Cubing an ideal", "[CubeNucube]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(193));
    QuadraticNumber<ZZ> quad_number1 = QuadraticNumber<ZZ>(quad_order1);
    QuadraticNumber<ZZ> quad_number2 = QuadraticNumber<ZZ>(quad_order1);

    CubeNucube<ZZ> cube_nucube_object = CubeNucube<ZZ>();
    cube_nucube_object.set_RelativeGenerator(quad_number1);
    quad_order1.set_cube_nucube(cube_nucube_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    red_plain_real_object.set_RelativeGenerator(quad_number2);
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    QuadraticIdealBase<ZZ> quad_ideal_base_cube = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(2), ZZ(13), ZZ(-3));

    cube(quad_ideal_base_cube, quad_ideal_base1);

    REQUIRE(quad_ideal_base_cube.get_a() == 6);
    REQUIRE(quad_ideal_base_cube.get_b() == 11);
    REQUIRE(quad_ideal_base_cube.get_c() == -3);

}
#endif
