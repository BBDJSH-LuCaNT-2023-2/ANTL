/**
 * @file zz_p_expo_test.cpp
 * @author Michael Jacobson
 * @brief Test program for exponentiation routines using zz_p as base type.
 */

#include <NTL/lzz_p.h>
#include <ANTL/Exponentiation/ExponentiationBinary.hpp>
#include <ANTL/Exponentiation/ExponentiationNAF.hpp>
#include <ANTL/Exponentiation/ExponentiationL2R.hpp>
#include <ANTL/Exponentiation/ExponentiationWNAF.hpp>
#include <ANTL/Exponentiation/ExponentiationYao.hpp>



NTL_CLIENT
using namespace ANTL;
	
int main (int argc, char **argv)
{
  zz_p a,b_bin,b_naf, b_l2r, b_wnaf, b_yao;
  ZZ n;

  // use GF(1073741827) for these tests
  zz_p::init(1073741827);

  // set base to be a random value in GF(1073741827)
  do {
    random(a);
  } while (IsOne(a));


  // generate random exponent of size 512 bits
  RandomLen (n, 512);

  cout << "Using:" << endl;
  cout << " p = " << zz_p::modulus() << endl;
  cout << " a = " << a << endl;
  cout << " n = " << n << endl;

  // initialize exponentiation classes
  ExponentiationBinary<zz_p> ebin;
  ExponentiationNAF<zz_p> enaf;
  ExponentiationL2R<zz_p> el2r;
  ExponentiationWNAF<zz_p> ewnaf;
  ExponentiationYao<zz_p> eyao;
 
  // compute a^n with available methods
  ebin.power(b_bin,a,n);

  enaf.initialize(a,n);
  enaf.power(b_naf,a,n);

  el2r.initialize(a);
  el2r.power(b_l2r,a,n);

  ewnaf.initialize(a,n,5);
  ewnaf.power(b_wnaf,a,n);

  eyao.power(b_yao,a,n);
 
  // check and output results
  cout << "a^n (binary) = " << b_bin << endl;
  cout << "a^n (naf)    = " << b_naf << endl;
  cout << "a^n (l2r)    = " << b_l2r << endl;
  cout << "a^n (wnaf)   = " << b_wnaf << endl;
  cout << "a^n (byao)   = " << b_yao << endl;

  if ((b_bin == b_naf) && (b_naf == b_l2r) && (b_wnaf == b_bin) && (b_bin == b_yao))
    cout << "RESULTS MATCH!" << endl;
  else
  {
    cout << "ERROR:  RESULTS DO NOT MATCH!" << endl;
    exit(1);
  }
}
