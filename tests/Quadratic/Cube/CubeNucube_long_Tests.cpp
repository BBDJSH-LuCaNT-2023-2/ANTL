#ifndef CUBE_NUCUBE_LONG_TEST
#define CUBE_NUCUBE_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Cube/CubeNucube.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("CubeNucube<long>: Cubing an ideal", "[CubeNucube]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(13);

    CubeNucube<long> cube_nucube_object = CubeNucube<long>();
    quad_order1.set_cube_nucube(cube_nucube_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base_cube = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, 9, 1);

    quad_order1.get_cube_nucube()->cube(quad_ideal_base_cube, quad_ideal_base1);

    REQUIRE(quad_ideal_base_cube.get_a() == 4913);
    REQUIRE(quad_ideal_base_cube.get_b() == 6299);
    REQUIRE(quad_ideal_base_cube.get_c() == 2019);

}
#endif
