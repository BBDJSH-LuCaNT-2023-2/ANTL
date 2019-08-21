#ifndef GUARD_cnftests_cpp
#define GUARD_cnftests_cpp

#include <cmath>
#include <boost/math/tools/polynomial.hpp>
#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/RealCubicNumberField.hpp"

// For floating point arithmetic error tolerance
const double DOUBLE_TOLERANCE = 0.0000001;

using std::abs;
using namespace NTL;
// 010-TestCase.cpp

//using namespace NTL;

// Let Catch provide main():

#include "../catch.hpp"


  int cubic_type;
  double a,b,c,d;

TEST_CASE("Root-finding function SolveP3, monic, double") {

  double rootArray[3];
  SECTION("Integral Roots"){
    a = 1; b = -2; c = -1; d = 2;
    //polynomial<long> P1{{d,c,b,a}};


    cubic_type = SolveP3<double>(rootArray, b,c,d);

    REQUIRE ( abs(rootArray[0] - (-1.0) ) < DOUBLE_TOLERANCE );
    REQUIRE ( abs(rootArray[2] - 1.0) < DOUBLE_TOLERANCE );
    REQUIRE ( abs(rootArray[1] - 2.0) <DOUBLE_TOLERANCE );
  }

  SECTION("Complex Roots"){
    a = 1; b = 1; c = 1; d = 4;
    cubic_type = SolveP3<double>(rootArray, b,c,d);

    REQUIRE ( abs(rootArray[0] - (-1.7429592) ) < DOUBLE_TOLERANCE );
    REQUIRE ( abs(rootArray[1] - (0.3714796) ) < DOUBLE_TOLERANCE );
    REQUIRE ( abs(abs(rootArray[2]) - (1.4686560) ) < DOUBLE_TOLERANCE );
  }

  SECTION("Not Monic, real"){
    a = 2; b = 1; c = -5; d = -2;
    cubic_type = SolveP3<double>(rootArray, b/a,c/a,d/a);

    REQUIRE ( abs(rootArray[0] - (-1.6485352) ) < DOUBLE_TOLERANCE );
    REQUIRE ( abs(rootArray[1] - (1.5419361) ) < DOUBLE_TOLERANCE );
    REQUIRE ( abs(rootArray[2] - (-0.3934009) ) < DOUBLE_TOLERANCE );
  }
} // end test case

boost::math::tools::polynomial<long> poly1{{-1,-3,1,1}};
boost::math::tools::polynomial<boost::multiprecision::mpz_int> poly2{{-1,-3,1,1}};
TEST_CASE("Discriminant function"){


  SECTION("Binary Cubic Form, discriminant long and mpz_int"){
    REQUIRE(discriminant_bcf(poly1) == 148);
    REQUIRE(discriminant_bcf(poly2) == 148);
  }
}

TEST_CASE("Constructor of RealCubicNF, and accessors"){
  RealCubicNumberField<long, double> Rufio(poly1);
  boost::math::tools::polynomial<long> * Rpoly;
  Rpoly = Rufio.get_defining_polynomial();

  SECTION("Testing Global accessors"){
    REQUIRE(Rufio.get_discriminant() == 148);
    REQUIRE( Rpoly->degree() == 3 );
    REQUIRE( (*Rpoly)[0] == -1 );
    REQUIRE( (*Rpoly)[1] == -3 );
    REQUIRE( (*Rpoly)[2] == 1 );
    REQUIRE(Rufio.is_real() == true);
    REQUIRE(abs(Rufio.get_root(0)) > 1.0);

  }
}
/*
TEST_CASE( "Multiply matrices", "[single-file]" ) {
  reducedIndexForm[0] = to_ZZ(1); //a
  reducedIndexForm[1] = to_ZZ(1); //b
  reducedIndexForm[2] = to_ZZ(1); //c
  reducedIndexForm[3] = to_ZZ(3); //d


  NMatrix M1, M2, M3;
  LargeNumber x,y,z,d;
  x = 1;
  y = 0;
  z = 0;
  d = 1;
  M2 = arithmeticMatrix(x,y,z,d);

  //Check this matrix is the identity
  REQUIRE(M2.denom == 1);
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
        if (i != j){
          REQUIRE(M2.mat[i][j] == 0);
        }
        else{
          REQUIRE(M2.mat[i][j] == 1);
        }
    }
  }

  x = RandomBnd(ZZ(47));
  y = RandomBnd(ZZ(22));
  z = RandomBnd(ZZ(15));
  d = RandomBnd(ZZ(19))+ZZ(1);

  M1 = arithmeticMatrix(x,y,z,d);
  M3 = arithmeticMatrix(x,y,z,d);
  SECTION(
    "Multiply by identity on the right"){
    MatrixMultiplication(M1,M2);

    bool entriesSame = true;

    for (int i = 0; i < 3; ++i){
      for (int j = 0; j < 3; ++j){
          if (M1.mat[i][j] != M3.mat[i][j]){
            entriesSame = false;
            break;
          }
      }
    }
    REQUIRE( entriesSame );
    REQUIRE(M1.denom == M3.denom);
  }

  SECTION("Multiply by identity on the left"){
    NMatrix M3;
    //copy left matrix to M3
    for (int i = 0; i < 3; ++i){
      for (int j = 0; j < 3; ++j){
          M3.mat[i][j] = M1.mat[i][j];
      }
    }
    M3.denom = M1.denom;

    //Compute M2 * M1
    MatrixMultiplication(M2,M1);

    bool entriesSame = true;

    for (int i = 0; i < 3; ++i){
      if (!entriesSame){
        break;
      }
      for (int j = 0; j < 3; ++j){
          if (M1.mat[i][j] != M3.mat[i][j]){
            entriesSame = false;
            break;
          }
      }
    }
    REQUIRE( entriesSame );
    REQUIRE(M2.denom == M3.denom);
  }


  SECTION("Multiply by two fixed matrices:"){
    x = 17;
    y = 5;
    z = 6;
    d = 12;
    M1 = arithmeticMatrix(x,y,z,d);
    x = 10;
    y = 5;
    z = 14;
    d = 6;
    M2 = arithmeticMatrix(x,y,z,d);

    REQUIRE(M1.mat[0][0]== 17);
    REQUIRE(M1.mat[1][0]== 5);
    REQUIRE(M1.mat[2][0]== 6);
    REQUIRE(M1.mat[0][1]== -18);
    REQUIRE(M1.mat[1][1]== 6);
    REQUIRE(M1.mat[2][1]== 5);
    REQUIRE(M1.mat[0][2]== -33);
    REQUIRE(M1.mat[1][2]== -23);
    REQUIRE(M1.mat[2][2]== 11);
    REQUIRE(M1.denom == 12);

    MatrixMultiplication(M1, M2);

    REQUIRE(M1.mat[0][0]== -382);
    REQUIRE(M1.mat[1][0]== -242);
    REQUIRE(M1.mat[2][0]== 239);
    REQUIRE(M1.mat[0][1]== -717);
    REQUIRE(M1.mat[1][1]== -379);
    REQUIRE(M1.mat[2][1]== -242);
    REQUIRE(M1.mat[0][2]== 9);
    REQUIRE(M1.mat[1][2]== -475);
    REQUIRE(M1.mat[2][2]== -621);
    REQUIRE(M1.denom == 72);
  }
}//end TEST_CASE

TEST_CASE( "Invert Matrix") {
  reducedIndexForm[0] = to_ZZ(1); //a
  reducedIndexForm[1] = to_ZZ(1); //b
  reducedIndexForm[2] = to_ZZ(1); //c
  reducedIndexForm[3] = to_ZZ(3); //d

  NMatrix M1, M2;
  LargeNumber x,y,z,d;
  x = 17;
  y = 5;
  z = 6;
  d = 12;

  SECTION("Inverting twice should return the original"){

    M1 = arithmeticMatrix(x,y,z,d);

    M2 = Inverse(M1);
    REQUIRE(M2.mat[0][0]==  1086);
    REQUIRE(M2.mat[1][0]== -1158);
    REQUIRE(M2.mat[2][0]==   -66);
    REQUIRE(M2.mat[0][1]==   198);
    REQUIRE(M2.mat[1][1]==  2310);
    REQUIRE(M2.mat[2][1]== -1158);
    REQUIRE(M2.mat[0][2]==  3672);
    REQUIRE(M2.mat[1][2]==  1356);
    REQUIRE(M2.mat[2][2]==  1152);
    REQUIRE(M2.denom == 3457);

    M1 = Inverse(M2);
    REQUIRE(M1.mat[0][0]== 17);
    REQUIRE(M1.mat[1][0]== 5);
    REQUIRE(M1.mat[2][0]== 6);
    REQUIRE(M1.mat[0][1]== -18);
    REQUIRE(M1.mat[1][1]== 6);
    REQUIRE(M1.mat[2][1]== 5);
    REQUIRE(M1.mat[0][2]== -33);
    REQUIRE(M1.mat[1][2]== -23);
    REQUIRE(M1.mat[2][2]== 11);
    REQUIRE(M1.denom == 12);
  }

  SECTION("Inverting a matrix with a common factor"){
    x = 2;
    y = 2;
    z = 2;
    d = 6;

    M1 = arithmeticMatrix(x,y,z,d);
    REQUIRE(M1.mat[0][0]== 1);
    REQUIRE(M1.mat[1][0]== 1);
    REQUIRE(M1.mat[2][0]== 1);
    REQUIRE(M1.mat[0][1]== -3);
    REQUIRE(M1.mat[1][1]== -1);
    REQUIRE(M1.mat[2][1]== 1);
    REQUIRE(M1.mat[0][2]== -6);
    REQUIRE(M1.mat[1][2]== -4);
    REQUIRE(M1.mat[2][2]== 0);
    REQUIRE(M1.denom == 3);
  }
}
*/

#endif
