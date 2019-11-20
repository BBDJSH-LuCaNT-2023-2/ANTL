#ifndef GENERAL_FUNCTIONS_H
#define GENERAL_FUNCTIONS_H

#include <iostream>
#include <complex>
#include <cmath>
#include <mpfr.h>
#include <mpc.h>


#include "../../include/ANTL/Cubic/generalFunctions.hpp"
#include <boost/multiprecision/gmp.hpp>

using namespace boost::multiprecision;
using boost::multiprecision::mpf_float_100;

using std::complex;




const double PI = 3.141592653589793;


// Sokal's method for evaluating the dedekind eta functions
// N is the number of terms to calculate.
complex<double> dedekindEta(complex<double>& z, long N){

  if (N == 0){
    return complex<double>(0.0, 0.0);
  }
  if (N == 1){
    return complex<double>(1.0, 0.0);
  }

  complex<double> eta = complex<double>(1.0, 0.0);
  complex<double> imaginary = complex<double>(0.0, 1.0);


  complex<double>x_value = exp(2*PI*imaginary*z); // e^{2 pi i z}

  complex<double> x_power = x_value;
  complex<double> term = complex<double>(1.0, 0.0);

  for (int j = 1; j < N; j++){
    term = -term;
    term *= x_power;
    term /= (1.0-x_power);

    x_power *= x_value;           // on the jth iteration, updates x_value^j -> x_value^{j+1}

    eta += term;
  }//end for

  return eta;
}









#endif // guard
