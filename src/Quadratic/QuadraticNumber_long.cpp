/**
 * @file QuadraticNumber_long.cpp
 * @author Michael Jacobson
 * @remark Quadratic number function specializations (long base type).
 */

#include <QuadraticNumber.hpp>


template <>
void
QuadraticNumber<long>::normalize()
{
	if(d < 0){
	  a = -a;
	  b = -b;
	  d = -d;
	}
	ZZ g = GCD(GCD(a,b),d);
	if (g != 1) {
	  a /= g;
	  b /= g;
	  d /= g;
	}
      }
}



template <>
void
QuadraticNumber<long>::isUnit()
{
  return ::abs(getNorm()) == 1; 
}
