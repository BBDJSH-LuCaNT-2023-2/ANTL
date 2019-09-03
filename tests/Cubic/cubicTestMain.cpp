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

  boost::qvm::mat<mpf_float_100,3,3> Matty = {{1,2,3,4,5,6,7,8,9}};
  boost::qvm::mat<double,3,3> Matilda = {{0,0,1,0,1,0,1,0,0}};
  boost::qvm::inverse(Matilda);
  std::cout << Matty.a[1][1]<< endl;
  for (int i = 0; i < 3; i++){
    cout << Matilda.a[i][0] << " " << Matilda.a[i][1] << " " << Matilda.a[i][2] << endl;
  }
  std::cout << boost::qvm::A22(Matty)<< endl;


  // testing boost polynomial
  polynomial<mpz_int> test_poly{{3,1,1,1}};
  polynomial<long> test_poly2{{1,1,-2,-1}};
  cout << test_poly[0]<< test_poly[1] << test_poly[2]<< test_poly[3] << endl;
  polynomial<mpz_int> * testptr = &test_poly;

  std::cout << "Pointer practice1" << (*testptr)[0] << (*testptr)[1] << std::endl;

  //ComplexCubicNumberField<long,long> Joe; //= ComplexCubicNumberField<long,long>(test_poly);
  //polynomial<long> * pptr = CNF.get_defining_polynomial();
  //cout << (pptr)->degree()<< endl;

  ANTL::QQ<long> rholder(1,2);
  RealCubicNumberField<long,double> Jonah(test_poly2);
  ComplexCubicNumberField<long,double> Jim(test_poly);
  cout << "discriminant: " << Jim.get_discriminant() << endl;
  cout << "real discriminant: " << Jonah.get_discriminant() << endl;

  CubicOrderNF<long, double> Odie(test_poly);
  cout << "order disc:  "<< Odie.get_discriminant() << endl;


  const CubicOrderNF<long, double> * odie_point = &Odie;
  long a = 0;
  //cout << ANTL::mul(a, 4L,5L) << endl;
  long ellie_coeffs[3] = {3,6,9};
  long elias_coeffs[3] = {1,0,0};
  long elmira_coeffs[3] = {3,3,3};
  CubicElementNF<long, double> Ellie( odie_point,  ellie_coeffs, 3);
  CubicElementNF<long, double> Elias( odie_point,  elias_coeffs, 1);
  CubicElementNF<long, double> Elmira( odie_point,  elmira_coeffs, 1);
  cout << "Ellie reduced:" << Ellie.get_u() << " " << Ellie.get_x() << " " << Ellie.get_y() << endl;
  mul(Ellie, Elias, Elmira);
  cout << "Ellie Values:" << Ellie.get_u() << " " << Ellie.get_x() << " " << Ellie.get_y() << endl;
  Ellie.trace(rholder);
  cout << "Ellie trace: "<< rholder.getNumerator()<< endl; add(Ellie, Elias, 100L);
  Ellie.assign(Elias);

  cout << "see below" << endl;
  assign(Elmira, Ellie);
  std::cout << (Elmira.get_order()->is_equal( (*Ellie.get_order()) ) ) << std::endl;
  cout << "see above" << endl;
  cout << "Elmira Values:" << Elmira.get_u() << " " << Elmira.get_x() << " " << Elmira.get_y() << endl;

  std::complex<double> z1 = (1.0l, 1.0);
  std::complex<NTL::RR> z2 = (RR(1.0), RR(1.0));
  cout << dedekindEta(z1, 0) << endl;
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

  return 0;
}
