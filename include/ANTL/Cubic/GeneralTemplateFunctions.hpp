#ifndef GENERAL_TEMPLATE_FUNCTIONS_HPP
#define GENERAL_TEMPLATE_FUNCTIONS_HPP

#include <boost/math/bindings/rr.hpp>
#include <boost/math/tools/polynomial.hpp>
#include "../common.hpp"
#include "../Arithmetic/QQ.hpp"

using namespace NTL;

using namespace ANTL;
using boost::math::tools::polynomial;
using boost::math::ntl::atan;
using boost::math::ntl::cos;

const double eps=1e-25;
const double PI = 3.141592653589793;
const	double TwoPi = 6.28318530717958648;
const RR DOUBLE_TOL(1e-25);


//const quad_float eps=1e-23;
//const quad_float PI =    3.14159265358979323846264;
//const	quad_float TwoPi = 6.28318530717958647692528;



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
/*
//////////////////////////////////////////////////////////////
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] � i*x[2], return 1
template<typename PP>
int SolveP3(NTL::RR* x,PP a,PP b,PP c){

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
  		PP t=r/sqrt(q3);
  		if( t<-1) t=-1;
  		if( t> 1) t= 1;
          t=acos(t);
          a/=3; q=-2*sqrt(q);
          x[0]=q*cos(t/3)-a;
          x[1]=q*cos((t+TwoPi)/3)-a;
          x[2]=q*cos((t-TwoPi)/3)-a;
          return(3);
      } else {
          //A =-pow(fabs(r)+sqrt(r2-q3),1./3);
          A =-root3<PP>(fabs(r)+sqrt(r2-q3));
  		if( r<0 ) A=-A;
  		B = (A==0? 0 : q/A);

  		a/=3;
  		x[0] =conv<RR>((A+B)-a);
          x[1] =conv<RR>(-0.5*(A+B)-a);
          x[2] = conv<RR>(0.5*sqrt(3.0)*(A-B) );
  		if(fabs(x[2])<eps) { x[2]=x[1]; return(2); }
          return(1);
      }
};*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template<typename T, typename PT>
int cardano(polynomial<T> const &poly, PT * roots){
  if (poly.degree() != 3){
    return -1;
  }else{
    QQ<T> p = QQ<T>(poly[1]);
    QQ<T> q = QQ<T>(poly[0]);                 // set as d
    QQ<T> q2;
    QQ<T> p3;
    T temp;
    PT realtemp1,realtemp2;
    PT disc, disc_root;
    const PT athird(PT(1)/PT(3));

    // The code below transforms a generic cubic into one of the form x^3 + px+q
    QQ<T> intermediate = QQ<T>(poly[2]);        // set as b
    sqr(intermediate, intermediate);            // b^2
    mul(temp, T(3), poly[3]);                   // 3a
    div(intermediate, intermediate, temp);      // b^2/3a
    sub(p, p, intermediate);


    mul(temp, T(2), poly[2]);
    mul(temp, temp, poly[2]);
    mul(temp, temp, poly[2]);
    intermediate.assign(temp);
    mul(temp, poly[3], poly[3]);
    mul(temp, temp, T(27));
    div(intermediate, intermediate, temp);
    add(q, q, intermediate);                     // q = d + (2b^3/27a^2)

    mul(temp, poly[2], poly[1]);              // bc
    intermediate.assign(temp);
    mul(temp, T(3), poly[3]);                // 3a
    div(intermediate, intermediate, temp);     // bc/3a
    sub(q, q, intermediate);                  // q = d + (2b^3/27a^2) - bc/3a

    if (! (IsOne(poly[3]) )){
      div(q, q, poly[3]);
      div(p, p, poly[3]);
    }
    // at this point we will be working with the depressed monic cubic
    // x^3 + px + q, whose coefficients are rational numbers

    div(q, q, T(2));                         // obtain q/2
    div(p, p, T(3));                         // obtain p/3

    sqr(q2, q);                              // (q/2)^2
    cube(p3, p);                             // (p/3)^3

    div(realtemp1, to<PT>(q2.getNumerator()), to<PT>(q2.getDenominator()));       //real value of (q/2)^2
    div(realtemp2, to<PT>(p3.getN()), to<PT>(p3.getD()));                         //real value of (p/3)^3

    add(disc, realtemp1, realtemp2);                                              // (q/2)^2 + (p/3)^3
    //std::cout << "Discrim " << disc << std::endl;


    // make sure that the discriminant is not 0
    if (disc < DOUBLE_TOL && disc > -DOUBLE_TOL){
      std::cout << "repeated roots, try a new polynomial" << std::endl;
      return 0;
    }
    // Checks if we are in the totally real case
    else if (disc< DOUBLE_TOL && disc < -DOUBLE_TOL){

      PT asixth(1);
      div(asixth, asixth, PT(6));                           //unfortunate method of obtaining 1/6

      PT arg(2); div(arg, arg, PT(3));
      PT argshift; ComputePi(argshift);
      mul(argshift, argshift, arg);                          //2pi/3

      disc = -disc;

      SqrRoot(disc_root, disc);                             // sqrt(Delta) * i

      div(realtemp1, to<PT>(q.getN()), to<PT>(q.getD()));
      realtemp1 = -realtemp1;                               // realtemp1 = -q/2

      mul(realtemp2, realtemp1, realtemp1);                  // (q/2)^2
      add(realtemp2, realtemp2, disc);                       // Delta^2 + (q/2)^2


      //
      pow(realtemp2, realtemp2, asixth);
      mul (realtemp2, realtemp2, PT(2));                    // complex radius of roots


      div(arg, disc_root, realtemp1);                      // this should be tan(phi)

      if (arg < DOUBLE_TOL){
        arg = -arg;
        realtemp2 = -realtemp2;
      }
      //arg = atan(arg).value();
      atan_val(arg, arg);
      div(arg, arg, PT(3));
      //std::cout<< disc_root << " " << realtemp1 << std::endl;

      roots[0] = cos(arg);
      //std::cout << roots[0] << " " << realtemp2<< std::endl;
      mul(roots[0], roots[0], realtemp2);

      add(arg, arg, argshift);
      roots[1] = cos(arg);
      mul(roots[1], roots[1], realtemp2);

      add(arg, arg, argshift);
      roots[2] = cos(arg);
      mul(roots[2], roots[2], realtemp2);

      div(realtemp1, to<PT>(poly[2]), PT(3));
      div(realtemp1, realtemp1, to<PT>(poly[3]));
      sub(roots[0], roots[0] , realtemp1);
      sub(roots[1], roots[1] , realtemp1);
      sub(roots[2], roots[2] , realtemp1);

      return 1;
    }
    else{

      // This is where the complex case is handled
      PT omega_re(-1); div(omega_re, omega_re, to<PT>(2));
      PT omega_im(3);
      SqrRoot(omega_im, omega_im);

      div(omega_im, omega_im, PT(2));         /// omega_re - omega_im = w, (3rd root of unity)
                                              // oemga_re + omega_im = w^2

      SqrRoot(disc_root, disc);

      div(realtemp1, to<PT>(q.getN()), to<PT>(q.getD()));

      //std::cout << "p " << p.getNumerator()<< "/"<< p.getDenominator() << std::endl;
      //std::cout << "q " << q.getNumerator()<< "/"<< q.getDenominator() << std::endl;
      realtemp1 = -realtemp1;      // -q/2

      add(realtemp1, realtemp1, disc_root); // -q/2 + sqrt(disc)

      if (realtemp1 < DOUBLE_TOL){
        realtemp1 = -realtemp1;
        pow(realtemp1,realtemp1, athird);
        realtemp1 = -realtemp1;
      }else{
        pow(realtemp1,realtemp1, athird);          // v1, one of the solutions
      }

      // to get w1, such that w1 +v1 = root1, we find -p/3v1
      div(realtemp2, to<PT>(p.getN()), to<PT>(p.getD()) );
      realtemp2 = -realtemp2;               // this is -p/3
      div(realtemp2, realtemp2, realtemp1); // -p/3v1         // this is w1


      add(roots[0], realtemp2, realtemp1);                // this should be x1

      //std::cout << "root0 " << roots[0] << std::endl;

      mul(roots[1], realtemp1, omega_re);                   //re(v1, omega)

      mul(roots[2], realtemp1, omega_im);                   // im(v1, omega)

      mul(realtemp1, realtemp2, omega_re);                  // re(w1, omega^2)

      add(roots[1], realtemp1, roots[1]);

      omega_im = -omega_im;
      mul(realtemp1, omega_im, realtemp2);
      add(roots[2], realtemp1, roots[2]);


      // shift back to the original roots
      div(realtemp1, to<PT>(poly[2]), PT(3));
      div(realtemp1, realtemp1, to<PT>(poly[3]));
      sub(roots[0], roots[0] , realtemp1);
      sub(roots[1], roots[1] , realtemp1);

      return 2;
    }
  } // close the else


}

////////////////////////////////////////////////////////////////////////////////
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
