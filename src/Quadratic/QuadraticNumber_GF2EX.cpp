/**
 * @file QuadraticNumber_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Quadratic number function specializations (GF2EX base type).
 */

#include <QuadraticNumber.hpp>

      /**
       * @brief Inverts the QuadraticNumber
       */
      void
      invert<GF2EX> ()
      {
	// ((a + b rho) / d)^-1 = (ad + bd hx - bd rho) / (a^2 + b^2 fx + a b hx)
	::mul(newB,b,d);

	::mul(temp,newB,QO->getH());
	::mul(newA,a,d);
	::add(newA,newA,temp);

	::sqr(newD,a);
	::sqr(temp,b);
	::mul(temp,temp,QO->getF());
	::sub(newD,newD,temp);
	::mul(temp,a,b);
	::mul(temp,temp,QO->getH());

	a = newA;
	b = newB;
	d = newD;

	normalize ();
      }
