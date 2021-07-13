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

using namespace std;
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
  long dis, a,b,c,d;
  char waste;

  ofstream outfile;
  outfile.open("pari-formatted-polynomials.txt");
  outfile << "data = [\\ \n";
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
    std::cout << dis << "  " << line.substr(pos1, pos2) << std::endl;
    std::cout << a << " " << b << " " << c << " " << d << std::endl;

    outfile << "("<< a << ")*x^3 + (" << b<< ")*x^2 + (" << c << ")*x + (" <<d<< "),\\ "<<std::endl;
    // now a,b,c,d are the coefficients of the IBCF, we can now do the work:


  };
  outfile << "];";


  inFile.close();
  return 0;
}
