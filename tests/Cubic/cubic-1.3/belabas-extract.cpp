//#include <boost/math/tools/polynomial.hpp>

//#include "../../../include/ANTL/Arithmetic/QQ.hpp"
//#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"

#include<fstream>
#include<iostream>
#include <sstream>
#include <string>
#include <functional>
#include <iterator>
#include <unordered_map>



#include "../../../include/ANTL/Cubic/generalFunctions.hpp"
#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../../include/ANTL/Cubic/CubicOrder.hpp"
#include "../../../include/ANTL/Cubic/CubicElement.hpp"
#include "../../../include/ANTL/Cubic/CubicIdeal.hpp"

#include "../../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp"
#include "../../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"
#include <boost/math/bindings/rr.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfi.hpp>
using namespace NTL;
using namespace ANTL;
using boost::math::tools::polynomial;
using boost::math::ntl::atan;


int main(int argc,  char *argv[]){

  std::string line;
  std::string delimiter1 = "[";
  std::string delimiter2 = "]";
  int pos1, pos2;
  ifstream inFile;
  inFile.open(argv[1]);
  if (!inFile) {
    cerr << "Unable to open file datafile.txt";
    exit(1);   // call system to stop
  }
  std::cout << "we are in" << std::endl;
  NTL::ZZ dis, a,b,c,d;
  char waste;
  while (std::getline(inFile, line)){
    std::istringstream iss(line);
    iss >> dis;
    pos1 = line.find(delimiter1);
    pos2 = line.find(delimiter2);

    std::istringstream iss1(line.substr(pos1, pos2));
    iss1 >> waste; iss1 >> a;
    iss1 >> waste; iss1 >> b;
    iss1 >> waste; iss1 >> c;
    iss1 >> waste; iss1 >> d;
    std::cout << dis << "  " << line.substr(pos1, pos2) << endl;
    std::cout << a << " " << b << " " << c << " " << d << endl;

    polynomial<ZZ> const test_poly{{d,c,b,a }};
    std::cout << test_poly[3] << "x^3 + " << test_poly[2]<< "x^2 + " << test_poly[1] << "x + " << test_poly[0]<< std::endl;
    // now a,b,c,d are the coefficients of the IBCF, we can now do the work:


  };


  inFile.close();
  return 0;
}
