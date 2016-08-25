/**
 * @file QuadraticNumber_ZZ.cpp
 * @author Michael Jacobson
 * @remark Quadratic number function specializations (ZZ base type).
 */

#include <QuadraticNumber.hpp>


template <>
void
QuadraticNumber<ZZ>::normalize()
{
	if(d < 0){
	  negate(a,a);
	  negate(b,b);
	  negate(d,d);
	}
	ZZ g = GCD(GCD(a,b),d);
	if (g != 1) {
	  div(a,a,g);
	  div(b,b,g);
	  div(d,d,g);
	}
      }
}



template <>
void
QuadraticNumber<ZZ>::isUnit()
{
  return ::abs(getNorm()) == 1; 
}
