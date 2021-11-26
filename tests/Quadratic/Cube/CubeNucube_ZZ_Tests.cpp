#ifndef CUBE_NUCUBE_ZZ_TEST
#define CUBE_NUCUBE_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Cube/CubeNucube.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("CubeNucube<ZZ>: Cubing an ideal", "[CubeNucube]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    CubeNucube<ZZ> cube_nucube_object = CubeNucube<ZZ>();
    quad_order1.set_cube_nucube(cube_nucube_object);

    ReducePlainReal<ZZ> red_plain_real_object = ReducePlainReal<ZZ>();
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);

    QuadraticIdealBase<ZZ> quad_ideal_base_cube = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(-9), ZZ(1));

    cube(quad_ideal_base_cube, quad_ideal_base1);

    REQUIRE(quad_ideal_base_cube.get_a() == ZZ(1));
    REQUIRE(quad_ideal_base_cube.get_b() == ZZ(3));
    REQUIRE(quad_ideal_base_cube.get_c() == ZZ(-1));

}
#endif
