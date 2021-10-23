#ifndef GENERAL_TEMPLATE_FUNCTIONS_HPP
#define GENERAL_TEMPLATE_FUNCTIONS_HPP

#include <boost/math/bindings/rr.hpp>
#include <boost/math/tools/polynomial.hpp>
#include "../Arithmetic/QQ.hpp"
#include "../common.hpp"
#include <functional>
#include <cstdint>

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
Type calc_discriminant(polynomial<Type> const &poly){
  return (poly[2]*poly[2]*poly[1]*poly[1]) + Type(18)*poly[0]*poly[1]*poly[2]*poly[3] -
    Type(4)*poly[3]*poly[1]*poly[1]*poly[1] - Type(4)*poly[2]*poly[2]*poly[2]*poly[0] -
    Type(27)*poly[3]*poly[3]*poly[0]*poly[0];
}


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
    QQ<T> p = QQ<T>(poly[1]);                 // set as c
    QQ<T> q = QQ<T>(poly[0]);                 // set as d
    QQ<T> q2;
    QQ<T> p3;
    T temp;
    PT realtemp1,realtemp2;
    PT disc, disc_root;
    const PT athird(PT(1)/PT(3));
    // The code below transforms a generic cubic into one of the form x^3 + px+q
    // based on the writeup from https://brilliant.org/wiki/cardano-method/#_=_
    QQ<T> intermediate = QQ<T>(poly[2]);        // set as b
    sqr(intermediate, intermediate);            // b^2
    mul(temp, T(3), poly[3]);                   // 3a
    div(intermediate, intermediate, temp);      // b^2/3a
    sub(p, p, intermediate);                    //c - (b^2/3a)


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
    std::cout << "Depressed form: p,q " << p << "  " << q << endl;
    div(q, q, T(2));                         // obtain q/2
    div(p, p, T(3));                         // obtain p/3

    sqr(q2, q);                              // (q/2)^2
    cube(p3, p);                             // (p/3)^3

    div(realtemp1, to<PT>(q2.getNumerator()), to<PT>(q2.getDenominator()));       //real value of (q/2)^2
    div(realtemp2, to<PT>(p3.getN()), to<PT>(p3.getD()));                         //real value of (p/3)^3

    add(disc, realtemp1, realtemp2);                                              // (q/2)^2 + (p/3)^3


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


      pow(realtemp2, realtemp2, asixth);
      mul (realtemp2, realtemp2, PT(2));                    // complex radius of roots


      div(arg, disc_root, realtemp1);                      // this should be tan(phi)

      if (arg < DOUBLE_TOL){
        arg = -arg;
        realtemp2 = -realtemp2;
      }

      atan_val(arg, arg);
      div(arg, arg, PT(3));

      roots[0] = cos(arg);

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

      std::cout << "p " << p.getNumerator()<< "/"<< p.getDenominator() << std::endl;
      std::cout << "q " << q.getNumerator()<< "/"<< q.getDenominator() << std::endl;
      realtemp1 = -realtemp1;      // -q/2


      // special case when p = 0
      // see the description from brilliant.org for reference
      if( IsZero(p.getN()) ){
        mul(realtemp1, realtemp1, PT(2));
        std::cout << "realtemp1 " << realtemp1 << std::endl;
        if (realtemp1 < DOUBLE_TOL){
          realtemp1 = -realtemp1;
          pow(realtemp1,realtemp1, athird);
          realtemp1 = -realtemp1;
        }else{
          pow(realtemp1,realtemp1, athird);                   // S
        }
        roots[0] = realtemp1;
        mul(roots[1], realtemp1, omega_re);                   //re(S, omega)

        mul(roots[2], realtemp1, omega_im);                   // im(S, omega)
      }else{
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

        mul(roots[1], realtemp1, omega_re);                   //re(v1, omega)

        mul(roots[2], realtemp1, omega_im);                   // im(v1, omega)

        mul(realtemp1, realtemp2, omega_re);                  // re(w1, omega^2)

        add(roots[1], realtemp1, roots[1]);

        omega_im = -omega_im;
        mul(realtemp1, omega_im, realtemp2);
        add(roots[2], realtemp1, roots[2]);
      }

      // shift back to the original roots
      div(realtemp1, to<PT>(poly[2]), PT(3));
      div(realtemp1, realtemp1, to<PT>(poly[3]));
      sub(roots[0], roots[0] , realtemp1);
      sub(roots[1], roots[1] , realtemp1);

      std::cout << "Cardano roots: " << std::endl;
      std::cout << roots[0] <<std::endl;
      std::cout << roots[1] <<std::endl;
      std::cout << roots[2] <<std::endl;

      return 2;
    }
  }


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

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/**
* @params variable result will hold the output, evaluand is what to evaluate x at
* a,b,c,d are coefficients of the cubic ax^3 + bx^2 + cx + d
*/
template<typename Type>
void eval_cubic_mod_p(Type & result, const Type & evaluand, const Type &a, const Type &b, const Type & c, const Type & d, const Type & p){
  // horner's method for evaluating a cubic

  MulMod(result, a, evaluand, p);

  AddMod(result, result, b,p);
  MulMod(result, result, evaluand,p);
  AddMod(result, result, c,p);
  MulMod(result, result, evaluand,p);
  AddMod(result, result, d,p);
}






// based on code  obtained from  http://www.cplusplus.com/forum/general/85177/
// see also http://www.cplusplus.com/reference/cmath/atan/
// Only designed to work with RR and doubles
template<typename PT>
PT myAtan(PT argument, long terms)
{
	PT individual_term, sum, tempvar, temp_pi;
  individual_term = 0.0;
	sum = 0.0;

  ComputePi(temp_pi);
	// special cases
  abs(tempvar, argument-1.0);
	if( tempvar < to<PT>(DOUBLE_TOL) ) return (temp_pi/4.0);
  abs(tempvar, argument+1.0);
	if( tempvar < to<PT>(DOUBLE_TOL) ) return (-temp_pi/4.0);

	if(terms > 0)
	{
	    if( (argument < -1.0) || (argument > 1.0) )
	    {
	        if( argument > 1.0 )
                sum = temp_pi/2.0;
            else
                sum = -temp_pi/2.0;

            individual_term = -1.0/argument;
            for(int j=1; j<=terms; j++)
            {
                sum += individual_term;
                individual_term *= -1.0*(2.0*j-1)/((2.0*j+1)*argument*argument);
            }
	    }
	    else// -1 < x < 1
	    {
	        sum = 0.0;
            individual_term = argument;
            for(long j=1; j<=terms; j++)
            {
                sum += individual_term;
                individual_term *= -1.0*(2.0*j-1)*argument*argument/(2.0*j+1);
            }
	    }
	}

	return sum;
}

//long getPrecision(const double & x){return 64;}
//long getPrecision(const RR & x){return RR::precision();}


/**
* @brief Computes the value of the Dedekind Eta function. This one seems slightly off from
the Magma version
* @param[out] result is the value of the dedekind eta function at z
* @param[in] z a complex number on which to evalue eta(z)
* @param[in] N the number of terms to compute using the method of Sokal[2002]
* @param[in] result is a complex number to store the final value
*/
template<typename PT>
void DedekindEtaSokal(complex<PT>& result, const complex<PT>& z, const long & N){
    PT real_zero, real_num;
    real_zero = 0;
    real_num = 1;
    if (N < 0){result = complex<PT>(real_zero, real_zero);}
    if (N == 0){
      result = complex<PT>(real_zero, real_zero);
    }
    if (N == 1){
      result =  complex<PT>(real_num, real_zero);
    }

    complex<PT> eta = complex<PT>(real_num, real_zero);
    complex<PT> imaginary = complex<PT>(real_zero, real_num);

    ComputePi(real_num);
    mul(real_num, real_num, 2);

    result = exp(real_num*imaginary*z/to<PT>(24.0)); // e^{2 pi i z}

    complex<PT> x_value = exp(real_num*imaginary*z);
    complex<PT> x_power = x_value;
    complex<PT> term = complex<PT>(real_num, real_zero);

    for (int j = 1; j < N; j++){
      term = -term;
      term *= x_power;
      term /= (real_num-x_power);

      x_power *= x_value;           // on the jth iteration, updates x_value^j -> x_value^{j+1}

      eta += term;
    }
    result *= eta;
}

/**
* @brief Computes the value of the Dedekind Eta function. Uses the same formula
* described in Magma
* @param[out] result is the value of the dedekind eta function at z
* @param[in] z a complex number on which to evalue eta(z)
* @param[in] N the number of terms to compute
* @param[in] result is a complex number to store the final value
*/
template<typename PT>
void DedekindEta(complex<PT>& result, const complex<PT>& z, const long & N){
    PT real_zero, real_num, alt;
    real_zero = 0;
    real_num = 1;

    if (N == 0){
      result = complex<PT>(real_zero, real_zero);
    }
    if (N == 1){
      result =  complex<PT>(real_num, real_zero);
    }

    complex<PT> eta = complex<PT>(real_num, real_zero);
    complex<PT> imaginary = complex<PT>(real_zero, real_num);

    ComputePi(real_num);
    mul(real_num, real_num, 2);

    result = exp(real_num*imaginary*z/to<PT>(24.0)); // e^{2 pi i z}

    complex<PT> term = complex<PT>(real_num, real_zero);
    alt = -1;
    for (long j = 1; j < N; j++){
      term = alt;
      term *= exp(real_num*imaginary*z*to<PT>(j*(3*j-1)/2) )+exp(real_num*imaginary*z*to<PT>(j*(3*j+1)/2) );
      mul(alt, alt, -1);
      eta += term;
    }


    result *= eta;


}

// our hash function! Used in the hash table construction in BSGSVoronoi
namespace std {

  template<>
  struct hash<NTL::ZZ> {
    public:

      hash(NTL::ZZ m){
        modulus = m;
      }

      std::size_t operator()(const NTL::ZZ & a){


          std::size_t endvalue;
          rem(value_holder, a, modulus);

          conv(endvalue, value_holder);
          return endvalue;

      }



    private:
      ZZ modulus;
      ZZ value_holder;

  };
}

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
  void set_modulus(NTL::ZZ m){
    modulus = m;
  }

  ZZ modulus;
};

inline std::string zToString(const ZZ &z) {
    std::stringstream buffer;
    buffer << z;
    return buffer.str();
}
inline std::string zToString(const long &z) {
    std::stringstream buffer;
    buffer << z;
    return buffer.str();
}




#endif
