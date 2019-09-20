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

/*
template<typename Type>
void becomeCanonical(Type B[3][3]){

  Type g, r, s;          // store gcd g and integer multipliers r,s which
                               // satisfy r*m_32 + s*m_33 = g

  Type dummy, j,k,       //dummy is extraneous, j,k satisfy jr-ks = 1
              q,                // transitory value
              minorDet;         // determinant of minor matrix

  Type canonMatrix[3][3]; //variable for the canonical matrix


  updateLattice(B);             // function which removes any common factors
                                // in the denominator and numerators of entries

  // This if statement checks whether the properties of canonical basis are already
  //    satisfied
    if(
      ( B.coefficientMatrix[2][0] == 0 ) &&
      ( 0 <= B.coefficientMatrix[0][0] && B.coefficientMatrix[0][0] < B.mainDenominator) &&
      ( 0 <= B.coefficientMatrix[0][1] && B.coefficientMatrix[0][1] < B.mainDenominator) &&
      ( 0 <= B.coefficientMatrix[1][1] && B.coefficientMatrix[1][1] < B.coefficientMatrix[1][0])
    ) //close if clause
    {
        // do nothing if the condition is satisfied
    }
    else{
    //computes the determinant of the lower-right 2x2 minor
    minorDet = (B.coefficientMatrix[1][0] * B.coefficientMatrix[2][1]
                - B.coefficientMatrix[1][1] * B.coefficientMatrix[2][0]);

    //From the NTL Library, EEA
    XGCD(g, r, s ,B.coefficientMatrix[2][0], B.coefficientMatrix[2][1]);

    XGCD(dummy, j,k,r,-s);    // jr - ks = 1, dummy is always 1

    q = -(k*B.coefficientMatrix[2][0] + j*B.coefficientMatrix[2][1]);
    q = q/g;         // computes transitory value

    dummy = k+r*q;  //reusing dummy to save a couple ops


    canonMatrix[0][0] =  (dummy) * B.coefficientMatrix[0][0] + (j+q*s)*B.coefficientMatrix[0][1];
    canonMatrix[1][0] =  -minorDet/g;  // denoted  -epsilon/g in cubic book
    canonMatrix[2][0] =  0;
    canonMatrix[0][1] = r * B.coefficientMatrix[0][0] +  s*B.coefficientMatrix[0][1];
    canonMatrix[1][1] = r * B.coefficientMatrix[1][0] +  s*B.coefficientMatrix[1][1];
    canonMatrix[2][1] = g;

    // if -minorDet is negative, we need to flip the middle row to get e/g, where
    //  e is the absolute value of minorDet
    if (canonMatrix[1][0] < 0){
        canonMatrix[1][0] = -canonMatrix[1][0];
        canonMatrix[0][0] = - canonMatrix[0][0];
    }


    //Ensures that canonMatrix[1][1] (this is m_23 in the cubic book) satisfies
    // 0 <= m_23 < e/g )
    while ( (canonMatrix[1][1] >= canonMatrix[1][0]) || (canonMatrix[1][1] < 0) ){
        if (canonMatrix[1][1] < 0){
            canonMatrix[0][1] += canonMatrix[0][0];
            canonMatrix[1][1] += canonMatrix[1][0];
        }
        else if (canonMatrix[1][1] >= canonMatrix[1][0]){
            canonMatrix[0][1] -= canonMatrix[0][0];
            canonMatrix[1][1] -= canonMatrix[1][0];
        }
    }

    //ensures that m_12 is in the appropriate range
    while ( (canonMatrix[0][0] >= B.mainDenominator) || (canonMatrix[0][0] < 0) ){
        if (canonMatrix[0][0] < 0){
            canonMatrix[0][0] += B.mainDenominator;
        }
        else if (canonMatrix[0][0] >= B.mainDenominator){
            canonMatrix[0][0] -= B.mainDenominator;
        }
    }
    //ensures m_13 is in the appropriate range
    while ( (canonMatrix[0][1] >= B.mainDenominator) || (canonMatrix[0][1] < 0) ){
        if (canonMatrix[0][1] < 0){
            canonMatrix[0][1] += B.mainDenominator;
        }
        else if (canonMatrix[0][1] >= B.mainDenominator){
            canonMatrix[0][1] -= B.mainDenominator;
        }
    }

    for (int j = 0; j < 2 ; ++j){
      for (int i = 0; i <3; ++i){
        B.coefficientMatrix[i][j] = canonMatrix[i][j];
      }
    }
    updateLattice(B);
  } // close the main else clause

} //close function becomeCanonical
*/
#endif
