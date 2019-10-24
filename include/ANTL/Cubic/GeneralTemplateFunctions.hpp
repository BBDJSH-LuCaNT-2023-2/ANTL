#ifndef GENERAL_TEMPLATE_FUNCTIONS_HPP
#define GENERAL_TEMPLATE_FUNCTIONS_HPP

#include <boost/math/tools/polynomial.hpp>
#include "../Arithmetic/QQ.hpp"

using namespace ANTL;
using boost::math::tools::polynomial;



const double eps=1e-14;
const double PI = 3.141592653589793;
const	double TwoPi = 6.28318530717958648;
const double DOUBLE_TOL = 0.0000001;

template<typename PP>
static PP _root3 ( PP x );

template<typename PP>
PP root3 ( PP x );

template<typename PP>
int SolveP3(PP* x,PP a,PP b,PP c);

template<typename Type>
Type calc_discriminant(polynomial<Type> const &poly);


template<typename PP>
static PP _root3 ( PP x ){

  PP s = 1.;
  while ( x < 1. )
  {
      x *= 8.;
      s *= 0.5;
  }
  while ( x > 8. )
  {
      x *= 0.125;
      s *= 2.;
  }
  PP r = 1.5;
  r -= 1./3. * ( r - x / ( r * r ) );
  r -= 1./3. * ( r - x / ( r * r ) );
  r -= 1./3. * ( r - x / ( r * r ) );
  r -= 1./3. * ( r - x / ( r * r ) );
  r -= 1./3. * ( r - x / ( r * r ) );
  r -= 1./3. * ( r - x / ( r * r ) );
  return r * s;
};



template<typename PP>
PP root3 ( PP x ){
  if ( x > 0 ) return _root3<PP>( x ); else
  if ( x < 0 ) return-_root3<PP>(-x ); else
  return 0.;
};

//////////////////////////////////////////////////////////////
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
template<typename PP>
int SolveP3(PP* x,PP a,PP b,PP c){

    // Step 1: Convert to the canonical form  x^3 + q*x + r = 0
  	PP a2 = a*a;
    //std::cout << a2 << std::endl;
    PP q  = (a2 - 3*b)/9;
    //std::cout << q << std::endl;
  	PP r  = (a*(2*a2-9*b) + 27*c)/54;
  	// equation x^3 + q*x + r = 0


    PP r2 = r*r;
  	PP q3 = q*q*q;
  	PP A;
    PP B;
  	if (r2 <= (q3 + eps)) {//<<-- FIXED!
  		PP t=r/ANTL::sqrt(q3);
  		if( t<-1) t=-1;
  		if( t> 1) t= 1;
          t=acos(t);
          a/=3; q=-2*ANTL::sqrt(q);
          x[0]=q*cos(t/3)-a;
          x[1]=q*cos((t+TwoPi)/3)-a;
          x[2]=q*cos((t-TwoPi)/3)-a;
          return(3);
      } else {
          //A =-pow(fabs(r)+sqrt(r2-q3),1./3);
          A =-root3<PP>(fabs(r)+ANTL::sqrt(r2-q3));
  		if( r<0 ) A=-A;
  		B = (A==0? 0 : q/A);

  		a/=3;
  		x[0] =(A+B)-a;
          x[1] =-0.5*(A+B)-a;
          x[2] = 0.5*ANTL::sqrt(3.0)*(A-B);
  		if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
          return(1);
      }
};
////////////////////////////////////////////////////////////////////////////////

template<typename Type>
Type discriminant_bcf(polynomial<Type> const &poly){

// keep in mind that the boost polynomial class has index i corresponding to the ith power
  return poly[2]*poly[2]*poly[1]*poly[1]
  + Type(18)*poly[3]*poly[2]*poly[1]*poly[0]
  - Type(4)*poly[3]*poly[1]*poly[1]*poly[1]
  - Type(4)*poly[2]*poly[2]*poly[2]*poly[0]
  - Type(27)*poly[3]*poly[3]*poly[0]*poly[0] ;

}


#endif
