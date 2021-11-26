#ifndef CUBE_NUCUBE_LONG_TEST
#define CUBE_NUCUBE_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Cube/CubeNucube.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("CubeNucube<long>: Cubing an ideal", "[CubeNucube]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(long(13));

    CubeNucube<long> cube_nucube_object = CubeNucube<long>();
    quad_order1.set_cube_nucube(cube_nucube_object);

    ReducePlainReal<long> red_plain_real_object = ReducePlainReal<long>();
    quad_order1.set_red_best(red_plain_real_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);

    QuadraticIdealBase<long> quad_ideal_base_cube = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, -9, 1);

    cube(quad_ideal_base_cube, quad_ideal_base1);

    REQUIRE(quad_ideal_base_cube.get_a() == 1);
    REQUIRE(quad_ideal_base_cube.get_b() == 3);
    REQUIRE(quad_ideal_base_cube.get_c() == -1);

}
#endif
