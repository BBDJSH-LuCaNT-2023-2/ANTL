#ifndef GUARD_test_etafunction_cpp
#define GUARD_test_etafunction_cpp

#include <cmath>
#include <ctime>
#include<fstream>
#include<iostream>
#include <sstream>
#include <string>
#include <boost/math/tools/polynomial.hpp>
#include <ANTL/Cubic/GeneralTemplateFunctions.hpp>
#include "../catch.hpp"

// constants
//precision level requirement for double arithmetic
// Note that the input precision
const double DOUBLE_TOLERANCE = 0.000001;


TEST_CASE("Eta Function for doubles"){
  std::complex<double> z;
  std::complex<double> result;
  std::complex<double> im(0.0,1.0);
  std::string myfile = "TestInput/test-eta-input.txt";
  SECTION("equal tests"){

    REQUIRE(im.imag() > 0 );
    std::string line;
    double real_in, imag_in, real_out, imag_out;

    ifstream inFile;
    inFile.open(myfile);
    if (!inFile) {
      cerr << "Unable to open file";
      exit(1);   // call system to stop
    }

    for(int i =0; i < 50; ++i){
      std::getline(inFile, line);
      std::istringstream iss(line);
      iss >>  real_in >> imag_in >> real_out >> imag_out;

      z = real_in + im*imag_in;
      DedekindEta(result, z, 100);
      REQUIRE(std::abs(result.imag()-imag_out) < DOUBLE_TOLERANCE);
      REQUIRE(std::abs(result.real()-real_out) < DOUBLE_TOLERANCE);
    }

    inFile.close();
  }
}
TEST_CASE("Eta Function for RR"){
  std::complex<RR> z;
  std::complex<RR> result;
  RR abdif;
  std::complex<RR> im(to<RR>(0.0),to<RR>(1.0));
  std::string myfile = "etavalues.txt";
  SECTION("equal tests"){

    REQUIRE(im.imag() > 0 );
    std::string line;
    RR real_in, imag_in, real_out, imag_out;

    ifstream inFile;
    inFile.open(myfile);
    if (!inFile) {
      cerr << "Unable to open file";
      exit(1);   // call system to stop
    }

    for(int i =0; i < 50; ++i){
      std::getline(inFile, line);
      std::istringstream iss(line);
      iss >>  real_in >> imag_in >> real_out >> imag_out;

      z = real_in + im*imag_in;
      DedekindEta(result, z, 100);
      abs(abdif, to<RR>(result.imag())-imag_out);
      REQUIRE(abdif < to<RR>(DOUBLE_TOLERANCE));

      abs(abdif, to<RR>(result.real())-real_out);
      REQUIRE(abdif < to<RR>(DOUBLE_TOLERANCE));
    }

    inFile.close();
  }
}


#endif
