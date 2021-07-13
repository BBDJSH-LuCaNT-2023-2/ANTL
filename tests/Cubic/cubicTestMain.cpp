/** \mainpage
* This is a main file used for testing purposes
*/

#include <iostream>
#include <functional>
#include <iterator>
#include <unordered_map>
// qvm matrix headers
//#include <boost/qvm/mat.hpp>
//#include <boost/qvm/mat_traits.hpp>
//#include <boost/qvm/mat_access.hpp>
//#include <boost/qvm/mat_operations.hpp>

#include <boost/multiprecision/mpfi.hpp>

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

#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <complex>
using namespace NTL;
using namespace ANTL;
using NTL::ZZ;
using NTL::RR;
using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using boost::multiprecision::mpfi_float;
using std::cout;
using std::endl;
NTL_CLIENT


/*
struct ZZEqual {
 bool operator()(const ZZ& lhs, const ZZ& rhs) const
 {
    return lhs == rhs;
 }
};

struct ZZHash {

  ZZHash(NTL::ZZ m = ZZ(100)){
    modulus = m;
  }

  std::size_t operator()(const ZZ& a) const {
    std::size_t endvalue;
    ZZ value_holder;
    NTL::rem(value_holder, a, modulus);

    conv(endvalue, value_holder);
    return endvalue;
  }


  ZZ modulus;
};
*/

int main(){
  //std::unique_ptr< MultiplyStrategyWilliams<long,double> > pointo;


  //std::shared_ptr<MultiplyStrategyWilliams<long, double>> ims = std::make_shared<MultiplyStrategyWilliams<long, double>>();


  //mpfi_float mpfi_a = 2;
  //mpfi_float::default_precision(1000);
  //std::cout << mpfi_float::default_precision() << std::endl;
  //std::cout << sqrt(mpfi_a) << std::endl; // print root-2
  //long a,b,r,s,g;
  //a = 3; b= 5;
  //XGCD(g, r,s,a,b);
  //cout << "XGCD test: "<< g << " = " << a << "*" << r <<  " + " << b << "*" << s << endl;
  //ANTL::QQ<ZZ> rholder(to_ZZ(1),to_ZZ(2));

/* COMPLEX LONG INT
  polynomial<mpz_int> test_poly{{2,1,1,1}};
  cout << " Poly (C): " << test_poly[0]<< " + " << test_poly[1] << "x + " << test_poly[2]<< " x^2 + " << test_poly[3] << "x^3" << endl;
  CubicOrder<long, double> * co_point; co_point = CubicOrder<long, double>::make_order(test_poly);
  CubicOrder<long, double> * Odie = co_point;

  cout << "rho1:  " << Odie->get_rho1() << "     rho2:  " << Odie->get_rho2() << endl;
  cout << " root1 " << Odie->get_root1() << " root2 " << Odie->get_root2() << " root3 " << Odie->get_root3() << std::endl;
  cout << "order disc:  "<< Odie->get_discriminant() << endl;
  cout << "----------------------------------------" << endl;

  double val;
  Odie->get_FundamentalUnits(0)->get_real_value(val);
  std::cout << "Fundamental Units " << std::endl;
  std::cout << Odie->get_FundamentalUnits(0)->get_u() << " " << Odie->get_FundamentalUnits(0)->get_x() << " " << Odie->get_FundamentalUnits(0)->get_y() << " Reg: "<< ANTL::log(val) << std::endl;
*/


  ZZ ibcf[4];
  /*

  ibcf[3] = 1;
  ibcf[2] = 1;
  ibcf[1] = -3;
  ibcf[0] = -1;

  ibcf[3] = 2;
  ibcf[2] = 3;
  ibcf[1] = 5;
  ibcf[0] = 9;

  ibcf[3] = 1;
  ibcf[2] = 4;
  ibcf[1] = -2;
  ibcf[0] = -2;

*/
  ibcf[3] = 1;
  ibcf[2] = -1;
  ibcf[1] = 0;
  ibcf[0] = 8;


  polynomial<ZZ> const test_poly{{ibcf[0],ibcf[1],ibcf[2],ibcf[3] }};
  cout << "test_poly: "  << test_poly << endl;

  std::cout << test_poly[3] << "x^3 + " << test_poly[2]<< "x^2 + " << test_poly[1] << "x + " << test_poly[0]<< std::endl;
  CubicOrder<ZZ, RR> * co_point; co_point = CubicOrder<ZZ, RR>::make_order(test_poly);
  CubicOrder<ZZ, RR> * Odie = co_point;
  //Odie->roots_swap_position(0,1);Odie->roots_swap_position(1,2);
  ZZ E,G;
  Odie->standard_form(E,G);
  polynomial<ZZ> const standard_poly{{G,E,ZZ(0),ZZ(1) }};
  cout << E << " " << G << endl;
  cout << calc_discriminant(standard_poly) << endl;

  //Odie->set_unit_strategy("BSGS");
  cout << "rho1:  " << Odie->get_rho1() << "     rho2:  " << Odie->get_rho2() << endl;
  cout << " root1 " << Odie->get_root1() << " root2 " << Odie->get_root2() << " root3 " << Odie->get_root3() << std::endl;
  cout << "order disc:  "<< Odie->get_discriminant() << endl;
  cout << "----------------------------------------" << endl;

  RR val;
  Odie->get_fundamental_unit(0)->get_real_value(val);

  std::cout << test_poly[3] << "*x^3 + " << test_poly[2]<< "*x^2 + " << test_poly[1] << "*x + " << test_poly[0]<< std::endl;

  std::cout << "Fundamental Units " << std::endl;
  std::cout << Odie->get_fundamental_unit(0)->get_u() << " " << Odie->get_fundamental_unit(0)->get_x() << " " << Odie->get_fundamental_unit(0)->get_y() \
  << " Reg: "<< Odie->get_regulator() << std::endl;





//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// Real Order Testing
/*
ZZ real_ibcf[4];
real_ibcf[3] = 1;
real_ibcf[2] = 4;
real_ibcf[1] = -3;
real_ibcf[0] = -3;

polynomial<ZZ> const real_poly{{real_ibcf[0],real_ibcf[1],real_ibcf[2],real_ibcf[3] }};

  cout << " Poly (R): " << real_poly[0]<< " + " << real_poly[1] << "x + " << real_poly[2]<< " x^2 + " << real_poly[3] << "x^3" << endl;
  CubicOrder<ZZ, RR> * ro_point; ro_point = CubicOrder<ZZ, RR>::make_order(real_poly);
  CubicOrder<ZZ, RR> * Odessa = ro_point;

  Odessa->roots_swap_position(0,1);Odessa->roots_swap_position(1,2);
  cout << "rho1:  " << Odessa->get_rho1() << "     rho2:  " << Odessa->get_rho2() << endl;
  cout << " root1 " << Odessa->get_root1() << " root2 " << Odessa->get_root2() << " root3 " << Odessa->get_root3() << std::endl;
  cout << "order disc:  "<< Odessa->get_discriminant() << endl;
  cout << "----------------------------------------" << endl;
    std::cout << "Fundamental Units: " << std::endl;
  std::cout << Odessa->get_fundamental_unit(0)->get_u() << " " << Odessa->get_fundamental_unit(0)->get_x() << " " << Odessa->get_fundamental_unit(0)->get_y() << std::endl;
  std::cout << Odessa->get_fundamental_unit(1)->get_u() << " " << Odessa->get_fundamental_unit(1)->get_x() << " " << Odessa->get_fundamental_unit(1)->get_y() << std::endl;
  std::cout << Odessa->get_regulator() << endl;
  //std::cout << "Pointer practice1" << (*testptr)[0] << (*testptr)[1] << std::endl;
  //polynomial<long> * pptr = CNF.get_defining_polynomial();
  //cout << (pptr)->degree()<< endl;
*/

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  //polynomial<mpz_int> * testptr = &test_poly;



  //RealCubicNumberField<long,double> Jonah(real_poly);
  //ComplexCubicNumberField<long,double> Jim(test_poly);
  //CubicOrderReal<long, double> * Odie = dynamic_cast<CubicOrderReal<long, double>*>(cuboy);


  //cout << "rho1:  " << Odie->get_rho1() << "     rho2:  " << Odie->get_rho2() << endl;
  //cout << " root1 " << Odie->get_root1() << " root2 " << Odie->get_root2() << " root3 " << Odie->get_root3() << std::endl;
  //cout << "order disc:  "<< Odie->get_discriminant() << endl;
  //cout << "----------------------------------------" << endl;

  //CubicOrder<long, double> * odie_ptr = &Odie;


/*
  long number1[3] = {1,0,0};
  long number2[3] = {0,1,0};
  long number3[3] = {0,0,1};

  long n1[3] = {7,1,1};
  long n2[3] = {7,1,1};
  CubicElement<long, double> el1( Odie,  n1, 1);
  CubicElement<long, double> el2( Odie,  n2, 1);
  CubicElement<long, double> holder( Odie,  number1, 1);
  mul(holder,el1,el2 );

  long e1[3] = {1,0,0};
  long e2[3] = {0,1,0};
  long e3[3] = {0,0,1};

  CubicElement<long, double> Ellie( Odie,  number1, 1);
  CubicElement<long, double> Elias( Odie,  number2, 1);
  CubicElement<long, double> Elmira( Odie,  number3, 1);

  CubicElement<long, double> base1( Odie,  e1, 1);
  CubicElement<long, double> base2( Odie,  e2, 1);
  CubicElement<long, double> base3( Odie,  e3, 1);


  mul(holder, Ellie, base2);


  CubicIdeal<long, double> Ideal1(Odie, Ellie, Elias, Elmira);
  double punc_lat[2][2];
  */
  //Ideal1.make_prepared();
  //Ideal1.make_voronoi_basis();

  //std::cout << Odie->get_FundamentalUnits()->get_u() << " " << Odie->get_FundamentalUnits()->get_x() << " " << Odie->get_FundamentalUnits()->get_y() << std::endl;


/*
  CubicIdeal<long, double> Ideal2(Odie, base1, base2, base3);
  cout << "Ideal1 = " << Ideal1.get_coeff(0,0) << "  " << Ideal1.get_coeff(0,1) << "  " << Ideal1.get_coeff(0,2) \
  << "   Ideal2 = " << Ideal2.get_coeff(0,0) << "  " << Ideal2.get_coeff(0,1) << "  " << Ideal2.get_coeff(0,2) << endl;
  cout << "         " << Ideal1.get_coeff(1,0) << "  " << Ideal1.get_coeff(1,1) << "  " << Ideal1.get_coeff(1,2) \
  << "           " << Ideal2.get_coeff(1,0) << "  " << Ideal2.get_coeff(1,1) << "  " << Ideal2.get_coeff(1,2) << endl;
  cout << "         " << Ideal1.get_coeff(2,0) << "  " << Ideal1.get_coeff(2,1) << "  " << Ideal1.get_coeff(2,2) \
  << "            " << Ideal2.get_coeff(2,0) << "  " << Ideal2.get_coeff(2,1) << "  " << Ideal2.get_coeff(2,2) << endl;
  cout << " " << endl;
  cout << "         d="<<Ideal1.get_gen1()->get_denom() << "              d="<<Ideal2.get_denom() << endl;

  cout << "-------------------" << endl;
  cout << "-------------------" << endl;

  Odie->set_mul_strategy(ims);
  mul(Ideal1, Ideal1, Ideal2);

  cout << "Ideal product :" << Ideal1.get_coeff(0,0) << " " << Ideal1.get_coeff(0,1) << " " << Ideal1.get_coeff(0,2) << endl;
  cout << "              :" << Ideal1.get_coeff(1,0) << " " << Ideal1.get_coeff(1,1) << " " << Ideal1.get_coeff(1,2) << endl;
  cout << "              :" << Ideal1.get_coeff(2,0) << " " << Ideal1.get_coeff(2,1) << " " << Ideal1.get_coeff(2,2) << endl;
  cout << "        denom :"<< Ideal1.get_gen1()->get_denom() << endl;
*/

/*
  ZZHash h2(ZZ(100));
  std::unordered_multimap<ZZ, ZZ, ZZHash, ZZEqual> testmap{100, h2, ZZEqual()};

  testmap.reserve(100);
  //{{1, ZZ(10)}, {2, ZZ(12)}, {1, ZZ(15)}};
  cout << testmap.size() << endl;
  cout << testmap.count(ZZ(1)) << endl;
  //cout << testmap.count(5) << endl;
  cout << testmap.size() << endl;
  testmap.emplace(std::make_pair(ZZ(1), ZZ(12)));

  auto its = testmap.equal_range(ZZ(1));
  for (auto it = its.first; it != its.second; ++it) {
    cout << (it->second) << endl;
    cout << it->first << endl;
    cout << it->first << '\t' << it->second << endl;
  }
  //std::hash<ZZ> h1(ZZ(10));
  std::cout << h2(ZZ(15)) << std::endl;
*/

  return 0;
};
