#ifndef GENERAL_FUNCTIONS_H
#define GENERAL_FUNCTIONS_H

#include <complex>
#include <boost/math/bindings/rr.hpp>
#include <boost/multiprecision/gmp.hpp>
using namespace boost::multiprecision;
using boost::multiprecision::mpf_float;
using namespace std;



//#define	TwoPi  6.28318530717958648
//const double eps=1e-14;



// Sokal's method for evaluating the dedekind eta functions
// N is the number of terms to calculate.
complex<double> dedekindEta(complex<double>& z, long N);




#include "../../../src/Cubic/generalFunctions.cpp"
#endif // guard
