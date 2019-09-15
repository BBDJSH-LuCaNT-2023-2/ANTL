/** \mainpage
* This is a main file used for testing purposes
*/

#include <iostream>

// qvm matrix headers
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>


#include "../../include/ANTL/Cubic/generalFunctions.hpp"
#include "../../include/ANTL/Cubic/GeneralTemplateFunctions.hpp"
#include "../../include/ANTL/Cubic/CubicNumberField.hpp"
#include "../../include/ANTL/Cubic/RealCubicNumberField.hpp"
#include "../../include/ANTL/Cubic/ComplexCubicNumberField.hpp"
#include "../../include/ANTL/Cubic/CubicOrderNF.hpp"
#include "../../include/ANTL/Cubic/CubicElementNF.hpp"
#include "../../include/ANTL/Cubic/CubicIdealNF.hpp"

#include "../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp"
#include "../../include/ANTL/Cubic/Multiplication/MultiplyStrategyWilliams.hpp"

#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
using namespace NTL;
using namespace ANTL;
using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using std::cout;
using std::endl;
NTL_CLIENT

int main(){

  IdealMultiplicationStrategy<long, double> * ims = new MultiplyStrategyWilliams<long,double>();
  std::cout << ims->testmember << " also hi "<< std::endl;
  long a,b,r,s,g;
  a = 3;
  b= 5;
  XGCD(g, r,s,a,b);
  cout << g << " = " << a << "*" << r <<  " + " << b << "*" << s << endl;
  // testing boost polynomial
  polynomial<mpz_int> test_poly{{4,1,1,1}};
  polynomial<long> test_poly2{{1,1,-2,-1}};
  polynomial<mpz_int> * testptr = &test_poly;
  cout << test_poly[0]<< test_poly[1] << test_poly[2]<< test_poly[3] << endl;
  std::cout << "Pointer practice1" << (*testptr)[0] << (*testptr)[1] << std::endl;

  //ComplexCubicNumberField<long,long> Joe; //= ComplexCubicNumberField<long,long>(test_poly);
  //polynomial<long> * pptr = CNF.get_defining_polynomial();
  //cout << (pptr)->degree()<< endl;

  ANTL::QQ<long> rholder(1,2);
  RealCubicNumberField<long,double> Jonah(test_poly2);
  ComplexCubicNumberField<long,double> Jim(test_poly);
  cout << "Complex discriminant: " << Jim.get_discriminant() << endl;
  cout << "real discriminant: " << Jonah.get_discriminant() << endl;

  CubicOrderNF<long, double> Odie(test_poly);
  cout << "rho11  " << Odie.get_rho1() << endl;
  cout << "rho2  " << Odie.get_rho2() << endl;
  cout << "order disc:  "<< Odie.get_discriminant() << endl;

  const CubicOrderNF<long, double> * odie_ptr = &Odie;
  //cout << ANTL::mul(a, 4L,5L) << endl;
  long number1[3] = {3,6,9};
  long number2[3] = {1,0,0};
  long number3[3] = {3,3,3};

  long e1[3] = {1,0,0};
  long e2[3] = {0,1,0};
  long e3[3] = {0,0,1};
  CubicElementNF<long, double> base1( odie_ptr,  e1, 3);
  CubicElementNF<long, double> base2( odie_ptr,  e2, 3);
  CubicElementNF<long, double> base3( odie_ptr,  e3, 3);

  CubicElementNF<long, double> Ellie( odie_ptr,  number1, 3);
  CubicElementNF<long, double> Elias( odie_ptr,  number2, 1);
  CubicElementNF<long, double> Elmira( odie_ptr,  number3, 1);

  CubicIdealNF<long, double> Ideal1(odie_ptr, Ellie, Elias, Elmira);
  CubicIdealNF<long, double> Ideal2(odie_ptr, base1, base2, base3);

  Odie.set_mul_strategy(ims);
  Odie.mul(Ideal1, Ideal1, Ideal2);
  cout << "Ellie reduced:" << Ideal1.get_gen1()->get_u() << " " << Ideal1.get_gen2()->get_u() << " " << Ideal1.get_gen3()->get_u() << endl;
  cout << "Ellie reduced:" << Ideal1.get_gen1()->get_x() << " " << Ideal1.get_gen2()->get_x() << " " << Ideal1.get_gen3()->get_x() << endl;
  cout << "Ellie reduced:" << Ideal1.get_gen1()->get_y() << " " << Ideal1.get_gen2()->get_y() << " " << Ideal1.get_gen3()->get_y() << endl;
  cout << "Ideal denom "<<Ideal1.get_denom() << endl;
  cout << "Ellie reduced:" << Ellie.get_u() << " " << Ellie.get_x() << " " << Ellie.get_y() << endl;

  cout << "Ellie Values:" << Ellie.get_u() << " " << Ellie.get_x() << " " << Ellie.get_y() << endl;

  //CubicIdeal<long,double> Ideal1;


  //std::complex<double> z1 = (1.0l, 1.0);
  //std::complex<NTL::RR> z2 = (RR(1.0), RR(1.0));
  //cout << dedekindEta(z1, 0) << endl;
  //cout << dedekindEta(z1, 1) << endl;
  //cout << dedekindEta(z1, 2) << endl;
  //std::complex<RR> z3 = z2 + RR(1);
  cout << test_poly[3] << endl;



  mpf_float_100 root_list[3];
  mpf_float_100 a1 = test_poly[2]/test_poly[3];
  mpf_float_100 a2 = test_poly[1]/test_poly[3];
  mpf_float_100 a3 = test_poly[0]/test_poly[3];
  //std::cout << SolveP3<mpf_float_100>(root_list, a1,a2,a3) << std::endl;
  cout << Jonah.get_root(0) << " " << Jonah.get_root(1) << " " <<  Jonah.get_root(2) << endl;

  //delete ims;
  return 0;
}
