/**
 * @file common.cpp
 * @author Michael Jacobson
 * @brief implementation of non-templated functions defined in common.hpp
 */

#include <ANTL/common.hpp>

namespace ANTL {

  void DivRem(long & q, long & r, long a, long b)
  {
    q = (long)std::floor((double)a/(double)b);
    r = a - q*b;
  }

  //
  // Quadratic residuosity functions
  //

  /*
   * Function: Jacobi_base
   * Purpose: tests to see if a is a square mod p
   */
  template < class T, class U > inline long Jacobi_base (const T & a, const T & n, T & P1, T & P2, U & c, U & q) {
    long jac = 1;

    P1 = a;
    P2 = n;

    while (deg (P1) > 0) {
      c = rep (LeadCoeff (P1));
      if (deg (P2) & 1)
        jac *= ANTL::Jacobi(c, q);

      if (jac == 0)
        return jac;

      MakeMonic (P1);

      if ((deg (P1) & 1) && (deg (P2) & 1))
        jac *= ANTL::Jacobi(q - 1, q);


      swap (P1, P2);
      P1 %= P2;
    }

    if (deg (P2) & 1) {
      c = rep (LeadCoeff (P1));
      jac *= ANTL::Jacobi(c, q);
    }

  return jac;
  }

  long Jacobi_base (const long & a, const long & n) {
    long t,k,d;
    long nn,aa,temp;

    nn = n;
    aa = a;

    t = 1;

    while (aa != 0) {
      k = 0;
      while (!(aa & 1)) {
    aa >>= 1;
    ++k;
      }

      d = (nn & 7);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if ((aa & 3) == 3 && (d & 3) == 3) t = -t;
      temp = aa;
      aa = nn;
      nn = temp;
      aa %= nn;
    }

    if (nn == 1)
      return t;
    else
      return 0;
  }

  long Jacobi_base (const long long & a, const long long & n) {
    long t,k,d;
    long long nn,aa,temp;

    nn = n;
    aa = a;

    t = 1;

    while (aa != 0) {
      k = 0;
      while (!(aa & 1)) {
    aa >>= 1;
    ++k;
      }

      d = (nn & 7);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if ((aa & 3) == 3 && (d & 3) == 3) t = -t;
      temp = aa;
      aa = nn;
      nn = temp;
      aa %= nn;
    }

    if (nn == 1)
      return t;
    else
      return 0;
  }

  long Jacobi_base(const ZZ & a, const ZZ & n) {
    ZZ aa, nn;
    long t, k, d;

    nn = n;
    aa = a;
    t = 1;

    while (aa != 0) {
      k = MakeOdd(aa);
      d = trunc_long(nn, 3);
      if ((k & 1) && (d == 3 || d == 5)) t = -t;

      if (trunc_long(aa, 2) == 3 && (d & 3) == 3) t = -t;
      swap(aa, nn);
      rem(aa, aa, nn);
    }

    if (nn == 1)
      return t;
    else
      return 0;
  }


  long Jacobi(const long & a, const long & n) {
    long temp = a % n;
    if (temp < 0)  temp += n;
    return Jacobi_base(temp,n);
  }


  long Jacobi(const long long & a, const long long & n) {
    long long temp = a % n;
    if (temp < 0)  temp += n;
    return Jacobi_base(temp,n);
  }


  long Jacobi (const long long & a, const long & n) {
    long temp = (long) a % n;
    if (temp < 0)  temp += n;
    return Jacobi_base(temp,n);
  }


  long Jacobi(const ZZ & a, const ZZ & n) {
    ZZ temp;
    rem(temp,a,n);
    return Jacobi_base(temp,n);
  }


  long Jacobi (const ZZ & a, const long & n) {
    long temp = rem(a,n);
    return Jacobi_base(temp,n);
  }

  long Jacobi (const ZZ_pX & a, const ZZ_pX & n) {
    ZZ_pX temp;
    rem(temp,a,n);
    if (IsZero(temp))  return 0;

    ZZ_pX P1, P2;
    ZZ c,q;
    q = ZZ_p::modulus ();

    return Jacobi_base (a, n, P1, P2, c, q);
  }

  long Jacobi (const zz_pX & a, const zz_pX & n) {
    zz_pX temp;
    rem(temp,a,n);
    if (IsZero(temp))  return 0;

    zz_pX P1, P2;
    long c, q;
    q = (long) zz_p::modulus ();

    return Jacobi_base (a, n, P1, P2, c, q);
  }

  long Jacobi (const ZZ_pEX & a, const ZZ_pEX & n) {
    ZZ_pEX temp;
    rem(temp,a,n);
    if (IsZero(temp))  return 0;

    ZZ_pEX P1, P2;
    ZZ_pX c;
    ZZ_pX q = ZZ_pE::modulus ();

    return Jacobi_base (temp, n, P1, P2, c, q);
  }

  long Jacobi (const zz_pEX & a, const zz_pEX & n) {
    zz_pEX temp;
    rem(temp,a,n);
    if (IsZero(temp))  return 0;

    zz_pEX P1, P2;
    zz_pX c;
    zz_pX q = zz_pE::modulus ();

    return Jacobi_base (temp, n, P1, P2, c, q);
  }

  /*
   * Function: Kronecker
   * Purpose: This function calculates the kronecker symbol
   * Inputs: int Delta - the descriminante of the quadratic order
   *         int p - the prime divisor.
   * Output: long - the return value of the kronecker symbol.
   */
  long Kronecker(const long & a, const long & n) {
    long CHI;

    if (!(n & 1) && !(a & 1)) // a and n are both even
      CHI = 0;

    else {
      bool neg = false;
      long t = n;
      long temp;

      while (!(t & 1)) {
        t >>= 1;
        neg = !neg;
      }

      temp = a % t;
      if (temp < 0)  temp += t;
      CHI = Jacobi_base(temp,t);
      if (neg && ((a & 7) == 5))
      CHI = -CHI;
    }

    return CHI;
  }

  long Kronecker(const long long & a, const long long & n) {
    long CHI;

    if (!(n & 1) && !(a & 1)) // D and p are both even
      CHI = 0;

    else {
      bool neg = false;
      long long t = n;
      long long temp;

      while (!(t & 1)) {
        t >>= 1;
        neg = !neg;
      }

      temp = a % t;
      if (temp < 0)  temp += t;
      CHI = Jacobi_base(temp,t);
      if (neg && ((a & 7) == 5))
        CHI = -CHI;
    }

    return CHI;
  }

  long Kronecker(const long long & a, const long & n) {
    long CHI;

    if (!(n & 1) && !(a & 1)) // D and p are both even
      CHI = 0;

    else {
      bool neg = false;
      long t = n;
      long temp;

      while (!(t & 1)) {
        t >>= 1;
        neg = !neg;
      }

      temp = a % t;
      if (temp < 0)  temp += t;
      CHI = Jacobi_base(temp,t);
      if (neg && ((a & 7) == 5))
        CHI = -CHI;
    }

    return CHI;
  }

  long Kronecker(const ZZ & a, const ZZ & n) {
    long CHI;

    if (!IsOdd(n) && !IsOdd(a)) // D and p are both even
      CHI = 0;

    else {
      bool neg = false;
      ZZ t = n;
      ZZ temp;
      int D8 = trunc_long(a,3);
      if (a < 0)
        D8 = (8 - D8);

      while (!IsOdd(t)) {
        t >>= 1;
        neg = !neg;
      }

      rem(temp,a,t);
      CHI = Jacobi_base(temp,t);
      if (neg && (D8 == 5))
        CHI = -CHI;
    }

    return CHI;
  }


  long Kronecker(const ZZ & a, const long & n) {
    long CHI;

    if (!(n & 1) && !IsOdd(a)) // a and n are both even
      CHI = 0;

    else {
      bool neg = false;
      long t = n;
      long temp;
      int D8 = trunc_long(a,3);
      if (a < 0)
        D8 = (8 - D8);

      while (!(t & 1)) {
        t >>= 1;
        neg = !neg;
      }

      temp = rem(a,t);
      CHI = Jacobi_base(temp,t);
      if (neg && (D8 == 5))
        CHI = -CHI;
    }

    return CHI;
  }

  long Kronecker(const ZZ_pX & a, const ZZ_pX & n) {
    return Jacobi(a,n);
  }


  long Kronecker(const zz_pX & a, const zz_pX & n) {
    return Jacobi(a,n);
  }

  long Kronecker(const ZZ_pEX & a, const ZZ_pEX & n) {
    return Jacobi(a,n);
  }

  long Kronecker(const zz_pEX & a, const zz_pEX & n) {
    return Jacobi(a,n);
  }

//   long Kronecker(const GF2EX & h, const GF2EX & a, const GF2EX & n) {
//     return Jacobi(h,a,n);
//   }


//integer negation using mask m. m must be 0 or -1
// m=0: x->x
// m=-1: x->-x
//ported from liboptarith, written by Maxwell Sayles
int64_t negate_using_mask(const uint64_t m, const uint64_t x){
  assert(m == 0 || m == (uint64_t)(-1));
  return (x ^ m) - m;

}

} // ANTL
