/**
 * @file zz_p_expo_test.cpp
 * @author Michael Jacobson
 * @brief Test program for exponentiation routines using zz_p as base type.
 */

#include <NTL/lzz_p.h>
#include <ANTL/Exponentiation/ExponentiationBinary.hpp>
#include <ANTL/DExponentiation/DExponentiationIL.hpp>
#include <ANTL/DExponentiation/DExponentiationJSF.hpp>

NTL_CLIENT
using namespace ANTL;
	
int main (int argc, char **argv)
{
  zz_p a,b,ca_bin,cb_bin,c_bin, c_il, c_jsf;
  ZZ m,n;

  // use GF(1073741827) for these tests
  zz_p::init(1073741827);

  // set base to be a random value in GF(5)
  do {
    random(a);
  } while (IsOne(a));

  do {
    random(b);
  } while (IsOne(b));

  // generate random exponent of size 128 bits
  RandomLen (m, 512);
  RandomLen (n, 512);

  cout << "Using:" << endl;
  cout << " p = " << zz_p::modulus() << endl;
  cout << " a = " << a << endl;
  cout << " b = " << b << endl;
  cout << " m = " << m << endl;
  cout << " n = " << n << endl;

  // initialize exponentiation classes
  ExponentiationBinary<zz_p> ebin;
  DExponentiationIL<zz_p> deil;
  DExponentiationJSF<zz_p> djsf;

  // compute a^nb^m with available methods
  ebin.power(ca_bin,a,m);
  ebin.power(cb_bin,b,n);
  mul(c_bin,ca_bin,cb_bin);

  deil.initialize(a,b,m,n,5,5);
  deil.power(c_il,a,b,m,n);

  djsf.initialize(a,b,m,n);
  djsf.power(c_jsf,a,b,m,n);

  // check and output results
  cout << "a^mb^n (naive) = " << c_bin << endl;
  cout << "a^mb^n (interleaving)    = " << c_il << endl;
  cout << "a^mb^n (joint sparse form) = " << c_jsf << endl;

  if ((c_bin == c_il) and (c_il == c_jsf))
    cout << "RESULTS MATCH!" << endl;
  else
    cout << "ERROR:  RESULTS DO NOT MATCH!" << endl;
}
