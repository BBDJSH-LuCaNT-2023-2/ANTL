#ifndef CUBE_PLAIN_LONG_TEST
#define CUBE_PLAIN_LONG_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Cube/CubePlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("CubePlain<long>: Cubing two ideals", "[CubePlain]") {

    QuadraticOrder<long> quad_order1 = QuadraticOrder<long>(13);

    CubePlain<long> cube_plain_object = CubePlain<long>();
    quad_order1.set_cube_plain(cube_plain_object);

    QuadraticIdealBase<long> quad_ideal_base1 = QuadraticIdealBase<long>(quad_order1);
    QuadraticIdealBase<long> quad_ideal_base_cube = QuadraticIdealBase<long>(quad_order1);

    quad_ideal_base1.assign(17, 9, 1);

    quad_order1.get_cube_plain()->cube(quad_ideal_base_cube, quad_ideal_base1);

    REQUIRE(quad_ideal_base_cube.get_a() == 4913);
    REQUIRE(quad_ideal_base_cube.get_b() == 6299);
    REQUIRE(quad_ideal_base_cube.get_c() == 2019);

}
#endif
