/**
 * @file FPRepresentation.hpp
 * @author Reginald Lybbert
 * @brief reduced w-near (f,p) representations of ideals in real quadratic fields
 */

#ifndef ANTL_FP_REPRESENTATION_H
#define ANTL_FP_REPRESENTATION_H

#include <ANTL/common.hpp>
#include <ANTL/Arithmetic/QQ.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

namespace ANTL
{
  template < class T > class QuadraticIdealBase;
  template < class T > class QuadraticOrder;
  template < class T > class QuadraticNumber;


  // declare templated friend functions
  template < class T >
  void inv (FPRepresentation<T> & C, const FPRepresentation<T> &A);

  template < class T >
  void mul (FPRepresentation<T> &C, const FPRepresentation<T> &A, const FPRepresentation<T> &B);

  template < class T >
  void sqr (FPRepresentation<T> &C, const FPRepresentation<T> &A);

  template < class T >
  void cube (FPRepresentation<T> &C, const FPRepresentation<T> &A);

  template < class T >
  bool operator == (const FPRepresentation<T> &A, const FPRepresentation<T> &B);

  template < class T >
  bool operator != (const FPRepresentation<T> &A, const FPRepresentation<T> &B);

  template < class T >
  std::istream & operator >> (std::istream & in, FPRepresentation<T> &A);

  template < class T >
  std::ostream & operator << (std::ostream & out, const FPRepresentation<T> &A);


//Still need to write WNEAR, prove inversion, rewrite everything to match formats....
template < class T > class fpRepresentation
{
  protected:
    QuadraticIdealBase b;
    ZZ d; 
    long k;
    long w;
    long p;
    
    void remove(ZZ bigT,T C,long s){
       ZZ m = round(abs((bigT<<(p+3-s))/C));
       long t = floor(m/8d) + 1;
       d = ceil((1<<(p+3+t))*d/m);
       k = k - t;
    };

    void wnear();

public:

   //Still need to prove this works:
   friend void inv(FPRepresentation<T> & C, const FPRepresentation<T> &A){
      C.w = A.w;
      C.p = A.p;
      inv(C.b,A.b);
      C.d = floor(1<<(2*A.p + 1)/d);
      C.k = -A.k - 1;
   }
   friend void mul(FPRepresentation<T> & C, const FPRepresentation<T> & A, const FPRepresentation<T> &B){
       assert(A.w == B.w  && A.p == B.p);
       C.w = A.w;
       C.p = A.p;
       NUCOMP(C.b, A.b, B.b);
       if(A.d*B.d <= 1<<(2*A.p + 1)){
           C.d = (A.d*B.d)>>A.p;
           C.k = A.k + B.k;
       }else{
           C.d = (A.d*B.d)>>(A.p+1);
           C.k = A.k + B.k + 1; 
       }
       s = NumBits((C.b.B)/(C.b.Q))+(A.p+4);
       T = C.b.A<<s + C.b.B*floor(sqrt(D)<<s)
       C.remove(T,C.b.C,s);
       C.wnear();
   }

   friend void sqr(FPRepresentation<T> &C, const FPRepresentation<T> &A){
       C.w = A.w;
       C.p = A.p;
       NUDUPL(C.b, A.b);
       if(A.d*A.d <= 1<<(2*A.p + 1)){
           C.d = (A.d*A.d)>>A.p;
           C.k = 2*A.k;
       }else{
           C.d = (A.d*A.d)>>(A.p+1);
           C.k = 2*A.k + 1; 
       }
       s = NumBits((C.b.B)/(C.b.Q))+(A.p+4);
       T = C.b.A<<s + C.b.B*floor(sqrt(D)<<s)
       C.remove(T,C.b.C,s);
       C.wnear();
   };

   friend void cube(FPRepresentation<T> &C, const FPRepresentation<T> &A){
       C.w = A.w;
       C.p = A.p;
       NUCUBE(C.b, A.b);
       if(A.d*A.d*A.d <= 1<<(3*A.p + 1)){
           C.d = (A.d*A.d*A.d)>>(2*A.p);
           C.k = 3*A.k;
       }else if(A.d*A.d*A.d <= 1<<(3*A.p + 2)){
           C.d = (A.d*A.d*A.d)>>(2*A.p+1);
           C.k = 3*A.k + 1; 
       }else{
           C.d = (A.d*A.d*A.d)>>(2*A.p+2);
           C.k = 3*A.k + 2; 
       }
       s = NumBits((C.b.B)/(C.b.Q))+(A.p+4);
       T = C.b.A<<s + C.b.B*floor(sqrt(D)<<s)
       C.remove(T,C.b.C,s);
       C.wnear();
   };

   friend bool operator == (const FPRepresentation<T> &A, const FPRepresentation<T> &B){
      return (A.b == B.b && A.d == B.d && A.k == B.k && A.w == B.w && A.p == B.p);
   }

  friend bool operator != (const FPRepresentation<T> &A, const FPRepresentation<T> &B){
      return !(A.b == B.b && A.d == B.d && A.k == B.k && A.w == B.w && A.p == B.p);
   }
   
  friend std::istream & operator >> < T > (std::istream & in,
					   QuadraticIdealBase<T> &A);

  friend std::ostream & operator << < T > (std::ostream & out,
					   const QuadraticIdealBase<T> &A);
};


