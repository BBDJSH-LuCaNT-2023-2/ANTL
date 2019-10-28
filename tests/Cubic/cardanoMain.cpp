#include <boost/math/tools/polynomial.hpp>
#include <boost/math/bindings/rr.hpp>
#include "../../include/ANTL/Arithmetic/QQ.hpp"
#include "../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
using namespace NTL;
using namespace ANTL;
using boost::math::tools::polynomial;
using boost::math::ntl::atan;


int main(){

  ZZ a(2);
  ZZ b(3);
  ZZ c(-4);
  ZZ d(-2);
  polynomial<ZZ> mypol({d,c,b,a}) ;
  //std::cout << mypol[3] <<std::endl;

  RR::SetOutputPrecision(20);
  RR angle;
  angle = 0.333;

  //std::cout << atan(angle) << std::endl;
  RR rootvector[3];

  cardano(mypol, rootvector);
  std::cout << rootvector[0] << "  " << rootvector[1] << "  " << rootvector[2] << std::endl;
  return 0;
}
