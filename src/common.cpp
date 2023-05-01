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
template <class T>
T negate_using_mask(const uint64_t m, const T x){
  assert(m == 0 || m == (uint64_t)(-1));
  return (x ^ m) - m;

}
template<>
int64_t negate_using_mask<int64_t>(const uint64_t m, const int64_t x){
  assert(m == 0 || m == (uint64_t)(-1));
  return (x ^ m) - m;
}

int64_t sub_with_mask(uint64_t & m, const int64_t & a, const int64_t & b){
  int64_t r;
  #if defined(__x86_64)
  asm("subq %3, %0\n\t"
      "sbbq %1, %1\n\t"  // %1 is either 0 or -1
      : "=r"(r), "=&r"(m)
      : "0"(a), "r"(b)
      : "cc");

#else
  m = a < b ? -1 : 0;
  r = a - b;
#endif  
return r;
}

void cond_swap2_s64(int64_t & u1, int64_t & u2, int64_t & v1, int64_t & v2){
  uint64_t m;
  int64_t d2 = sub_with_mask(m, u2, v2);
  int64_t d1 = (u1 - v1) & m;
  d2 &= m;
  u1 -= d1;
  u2 -= d2;
  v1 += d1;
  v2 += d2;
}

uint64_t cond_swap3_s64(int64_t & u1,
				      int64_t & u2,
				      int64_t & u3,
				      int64_t & v1,
				      int64_t & v2,
				      int64_t & v3){
                  
  uint64_t m;
  int64_t d3 = sub_with_mask(m, u3, v3);
  int64_t d1 = (u1 - v1) & m;
  int64_t d2 = (u2 - v2) & m;
  d3 &= m;
  u1 -= d1;
  u2 -= d2;
  u3 -= d3;
  v1 += d1;
  v2 += d2;
  v3 += d3;
  return m;
}

int msb_u64(uint64_t x){
  #if defined(__x86_64)
  int64_t k = -1;
  asm("bsrq %1, %0\n\t"
      : "=r"(k)
      : "r"(x), "0"(k)
      : "cc");
  return k;
#else
  // a binary search approach to finding the most significant set bit
  int n = 0;
  if (x == 0) return -1;
  if (x > 0xFFFFFFFFULL) { n += 32; x >>= 32; }
  if (x > 0xFFFF) { n += 16; x >>= 16; }
  if (x > 0xFF) { n += 8;  x >>= 8; }
  if (x > 0xF) { n += 4;  x >>= 4; }
  if (x > 0x7) { n += 2;  x >>= 2; }
  if (x > 0x3) { n += 1;  x >>= 1; }
  if (x > 0x1) { n ++; }
  return n;
#endif
}

void ZZToMpz(const ZZ & A, mpz_t & a){
  //get array of mpz limbs from ZZ
  const mp_limb_t * zz_limbs = ZZ_limbs_get(A);
  //get number of limbs in ZZ 
  //the sign matters when reconstructing the mpz, see below
  long num_limbs = A.size(); 
  

  //get a pointer to a writable array of mp_limb_t s belonging to the mpz type
  mp_limb_t* mpz_limbs = mpz_limbs_write(a, num_limbs);
  for(int i=0; i<num_limbs;i++){
    mpz_limbs[i] = zz_limbs[i];
  }
  //mpz_limbs_finish uses the sign of num_limbs to set the sign on the mpz_t value  
  mpz_limbs_finish(a, sign(A)*num_limbs);
}

void MpzToZZ(const mpz_t & a, ZZ & A){
  //Get limbs from mpz_t
  const mp_limb_t * mpz_limbs = mpz_limbs_read(a);
  //Get number of limbs
  size_t num_limbs = mpz_size(a);
  //get mpz_t sign 
  int sign = mpz_sgn(a);

  //set the limbs of the ZZ object
  ZZ_limbs_set(A, mpz_limbs, num_limbs);
  //this seems slow
  A = A*ZZ(sign);

}
  
} // ANTL
