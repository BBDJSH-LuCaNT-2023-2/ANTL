#ifndef GUARD_test_voronoi_complex_cpp
#define GUARD_test_voronoi_complex_cpp

#include <cmath>
#include <ctime>
#include<fstream>
#include<iostream>
#include <sstream>
#include <string>
#include <boost/math/tools/polynomial.hpp>
#include <boost/math/bindings/rr.hpp>


#include <ANTL/Cubic/generalFunctions.hpp>
#include <ANTL/Cubic/GeneralTemplateFunctions.hpp>
#include <ANTL/Cubic/CubicNumberField.hpp>
#include <ANTL/Cubic/RealCubicNumberField.hpp>
#include <ANTL/Cubic/ComplexCubicNumberField.hpp>
#include <ANTL/Cubic/CubicOrder.hpp>
#include <ANTL/Cubic/CubicElement.hpp>
#include <ANTL/Cubic/CubicIdeal.hpp>
#include <ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp>
#include <ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp>
#include "../catch.hpp"

// constants
//precision level requirement for double arithmetic
// Note that the input precision
const double DOUBLE_TOLERANCE = 0.000001;


TEST_CASE("ZZ, RR Complex Voronoi tests"){
  std::string myfile = "TestInput/test-input-complex-regulator.txt";
  ZZ a,b,c,d, disc;
  RR test_regulator, temp;
  SECTION("Regulator Comparisons"){
    std::string line;
    double real_in, imag_in, real_out, imag_out;

    ifstream inFile;
    inFile.open(myfile);
    if (!inFile) {
      cerr << "Unable to open file";
      exit(1);   // call system to stop
    }

    for(int i =0; i < 100; ++i){
      std::getline(inFile, line);
      std::istringstream iss(line);
      iss >>  disc >> d >> c >> b >> a >> test_regulator;

      polynomial<ZZ> const real_poly{{d,c,b,a }};
      CubicOrder<ZZ, RR> * real_order_ptr = CubicOrder<ZZ,RR>::make_order(real_poly);
      abs(temp, real_order_ptr->get_regulator() - test_regulator);
      REQUIRE( temp < to<RR>(DOUBLE_TOLERANCE));

    }

    inFile.close();
  }
}



#endif
