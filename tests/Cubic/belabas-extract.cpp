//#include <boost/math/tools/polynomial.hpp>

//#include "../../../include/ANTL/Arithmetic/QQ.hpp"
//#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"

//#define DEBUG
//#define DEBUG2
//#define DEBUGVORONOI

#include<fstream>
#include<iostream>
#include <sstream>
#include <string>
#include <functional>
#include <iterator>
#include <unordered_map>
#include <vector>



#include "../../include/ANTL/Cubic/generalFunctions.hpp"
#include "../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../include/ANTL/Cubic/CubicOrder.hpp"
#include "../../include/ANTL/Cubic/CubicElement.hpp"
#include "../../include/ANTL/Cubic/CubicIdeal.hpp"

#include "../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp"
#include "../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"
#include <boost/math/bindings/rr.hpp>
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/multiprecision/mpfi.hpp>
using namespace NTL;
using namespace ANTL;
using NTL::ZZ;
using NTL::RR;
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

  ofstream outfile;
  outfile.open("bin/regulators.txt");

  NTL::ZZ dis, a,b,c,d;
  char waste;
  int num = 0;
  int start = 0;
  int end = 100;

  while (std::getline(inFile, line) && (num < end )){
    if(num >= start){
    std::istringstream iss(line);
    iss >> dis;
    pos1 = line.find(delimiter1);
    pos2 = line.find(delimiter2);

    std::istringstream iss1(line.substr(pos1, pos2));
    iss1 >> waste; iss1 >> a;
    iss1 >> waste; iss1 >> b;
    iss1 >> waste; iss1 >> c;
    iss1 >> waste; iss1 >> d;
      std::cout << "\nFIELD DATA: " << std::endl;
    std::cout << dis << "  " << line.substr(pos1, pos2) << endl;
    //std::cout << a << " " << b << " " << c << " " << d << endl;

    polynomial<ZZ> const real_poly{{d,c,b,a }};
    std::cout << "f(x) = "<<real_poly[3] << "x^3 + " << real_poly[2]<< "x^2 + " << real_poly[1] << "x + " << real_poly[0]<< std::endl;
    // now a,b,c,d are the coefficients of the IBCF, we can now do the work:

    //CubicOrderComplex<ZZ, RR> Order1(real_poly);
    CubicOrder<ZZ, RR> * ro_point = CubicOrder<ZZ,RR>::make_order(real_poly);

    CubicOrder<ZZ, RR> * Ord1 = ro_point;

    cout << "Disc of Order:  "<< Ord1->get_discriminant() << endl;
    //Ord1->roots_swap_position(1,2);//Ord1->roots_swap_position(1,2);
    std::cout << Ord1->get_root1() << " " << Ord1->get_root2() << " " << Ord1->get_root3() << std::endl;
    //cout << "conjugate elements" << Ord1->get_conjugate_bases(0,0) << " " << Ord1->get_conjugate_bases(0,1)<< endl;
    //cout << "conjugate elements" << Ord1->get_conjugate_bases(1,0) << " " << Ord1->get_conjugate_bases(1,1)<< endl;
    cout << "               " << endl;

    std::cout << "Calculating Fund. Units: " << std::endl;
    std::cout << "[" <<Ord1->get_fundamental_unit(0)->get_u() << ", " << Ord1->get_fundamental_unit(0)->get_x() << " " \
    << Ord1->get_fundamental_unit(0)->get_y() << "] ,  [" << Ord1->get_fundamental_unit(1)->get_u() << " "\
    << Ord1->get_fundamental_unit(1)->get_x() << " " << Ord1->get_fundamental_unit(1)->get_y() << "]"<<std::endl;
    outfile<< "Disc " << dis << " : FieldPoly: " << line.substr(pos1, pos2) << "  :  REGULATOR: "<< Ord1->get_regulator() << std::endl;
    //std::cout << "Pointer practice1" << (*testptr)[0] << (*testptr)[1] << std::endl;

    cout << "                                        " << endl;
    cout << "                                        " << endl;
    }
    num++;
  };

  outfile.close();
  inFile.close();
  return 0;
}
