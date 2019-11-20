#include <boost/math/tools/polynomial.hpp>
#include <boost/math/bindings/rr.hpp>
#include "../../../include/ANTL/Arithmetic/QQ.hpp"
#include "../../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"

#include<fstream>
#include<iostream>
#include <sstream>
#include <string>
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
  NTL::ZZ a;
  while (std::getline(inFile, line)){
    std::istringstream iss(line);
    iss >> a;
    pos1 = line.find(delimiter1);
    pos2 = line.find(delimiter2);

    std::cout << a << "  " << line.substr(pos1, pos2) << endl;
  };


  inFile.close();
  return 0;
}
