#ifndef CUBE_PLAIN_ZZ_TEST
#define CUBE_PLAIN_ZZ_TEST

#include "../../catch.hpp"
#include <ANTL/Quadratic/Cube/CubePlain.hpp>

using namespace NTL;
using namespace ANTL;

TEST_CASE("CubePlain<ZZ>: Cubing two ideals", "[CubePlain]") {

    QuadraticOrder<ZZ> quad_order1 = QuadraticOrder<ZZ>(ZZ(13));

    CubePlain<ZZ> cube_plain_object = CubePlain<ZZ>();
    quad_order1.set_cube_plain(cube_plain_object);

    QuadraticIdealBase<ZZ> quad_ideal_base1 = QuadraticIdealBase<ZZ>(quad_order1);
    QuadraticIdealBase<ZZ> quad_ideal_base_cube = QuadraticIdealBase<ZZ>(quad_order1);

    quad_ideal_base1.assign(ZZ(17), ZZ(9), ZZ(1));

    quad_order1.get_cube_plain()->cube(quad_ideal_base_cube, quad_ideal_base1);

    REQUIRE(quad_ideal_base_cube.get_a() == ZZ(4913));
    REQUIRE(quad_ideal_base_cube.get_b() == ZZ(6299));
    REQUIRE(quad_ideal_base_cube.get_c() == ZZ(2019));

}
#endif
