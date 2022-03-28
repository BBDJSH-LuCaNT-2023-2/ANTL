/**
 * @file xgcd_plain.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of XGCD_LEFT_PLAIN, XGCD_PARTIAL_PLAIN, and
 * XGCD_PARTIAL_REDUCE_PLAIN
 */

#include <ANTL/XGCD/xgcd_plain.hpp>

//
// XGCD_PLAIN
//

template <>
void XGCD_PLAIN(GF2EX & G, GF2EX & X, GF2EX & Y, const GF2EX & A, const GF2EX & B)
{
  XGCD_PLAIN_work<GF2E,GF2EX>(G,X,Y,A,B);
}

template <>
void XGCD_PLAIN(zz_pX & G, zz_pX & X, zz_pX & Y, const zz_pX & A, const zz_pX & B)
{
  XGCD_PLAIN_work<zz_p,zz_pX>(G,X,Y,A,B);
}

template <>
void XGCD_PLAIN(ZZ_pX & G, ZZ_pX & X, ZZ_pX & Y, const ZZ_pX & A, const ZZ_pX & B)
{
  XGCD_PLAIN_work<ZZ_p,ZZ_pX>(G,X,Y,A,B);
}

template <>
void XGCD_PLAIN(zz_pEX & G, zz_pEX & X, zz_pEX & Y, const zz_pEX & A, const zz_pEX & B)
{
  XGCD_PLAIN_work<zz_pE,zz_pEX>(G,X,Y,A,B);
}

template <>
void XGCD_PLAIN(ZZ_pEX & G, ZZ_pEX & X, ZZ_pEX & Y, const ZZ_pEX & A, const ZZ_pEX & B)
{
  XGCD_PLAIN_work<ZZ_pE,ZZ_pEX>(G,X,Y,A,B);
}



//
// XGCD_LEFT_PLAIN
//

template <>
void XGCD_LEFT_PLAIN(GF2EX & G, GF2EX & X, const GF2EX & A, const GF2EX & B)
{
  XGCD_LEFT_PLAIN_work<GF2E,GF2EX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PLAIN(zz_pX & G, zz_pX & X, const zz_pX & A, const zz_pX & B)
{
  XGCD_LEFT_PLAIN_work<zz_p,zz_pX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PLAIN(ZZ_pX & G, ZZ_pX & X, const ZZ_pX & A, const ZZ_pX & B)
{
  XGCD_LEFT_PLAIN_work<ZZ_p,ZZ_pX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PLAIN(zz_pEX & G, zz_pEX & X, const zz_pEX & A, const zz_pEX & B)
{
  XGCD_LEFT_PLAIN_work<zz_pE,zz_pEX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PLAIN(ZZ_pEX & G, ZZ_pEX & X, const ZZ_pEX & A, const ZZ_pEX & B)
{
  XGCD_LEFT_PLAIN_work<ZZ_pE,ZZ_pEX>(G,X,A,B);
}



//
// Partial Euclidean algorithm (for NUCOMP and fast reduce)
//


// no flag for GF2EX version
//  - Assumes that R1 is reduced modulo R2

void
XGCD_PARTIAL_PLAIN(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound)
{
#ifdef TRACE_XGCD
  cout << "--> IN XGCD_PARTIAL_PLAIN" << endl;
#endif
  long e = deg(R2) + 1;
  GF2EX q(INIT_SIZE, e), r(INIT_SIZE,e);

  clear(C2);
  set(C1);

  while (deg (R1) > bound)
    {
      DivRem (q, R2, R2, R1);
      swap(R2,R1);

      // r = C2 + q C1
      mul(r,q,C1);
      add(r,C2,r);
      C2 = C1;
      C1 = r;
    }
}



void
XGCD_PARTIAL_REDUCE_PLAIN(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even)
{
#ifdef TRACE_XGCD
  cout << "--> IN XGCD_PARTIAL_REDUCE_PLAIN" << endl;
#endif
  long e = deg(R2) + 1;
  GF2EX q(INIT_SIZE, e), r(INIT_SIZE,e);

  set(B2);
  clear(B1);

  while (deg (R1) > bound || (deg(R1) == bound && !even))
    {
      DivRem (q, R2, R2, R1);
      swap(R2,R1);

      // r = B2 + q B1
      mul(r,q,B1);
      add(r,B2,r);
      B2 = B1;
      B1 = r;
    }
}


//Plain XGCD from liboptarith 
//Written by Maxwell Sayles 
void 
XGCD_PLAIN(int64_t & g, int64_t & x, int64_t & y, const int64_t & a, const int64_t & b){

  int64_t t_a = a;
  int64_t t_b = b;

  uint64_t sm = a >> 63;
  uint64_t sn = b >> 63;
  t_a = ANTL::negate_using_mask(sm, a);
  t_b = ANTL::negate_using_mask(sn, b);
    
  int64_t m = 0;
  int64_t n = 1;
  int64_t u = 1;
  int64_t v = 0;
  
  if (t_a == 0) {
    x = 1;
    y = 0;
    g = b;
    return;
  }
  if (t_b == 0) {
    x = 0;
    y = 1;
    g = a;
    return;
  }

#if defined(__x86_64)
  asm("0:\n\t"
      "movq %0, %%rax\n\t"
      "xorq %%rdx, %%rdx\n\t"
      "divq %1\n\t"
      
      "movq %1, %0\n\t"
      "movq %%rdx, %1\n\t"
      
      "movq %%rax, %%rdx\n\t"
      "imul %4, %%rax\n\t"
      "imul %5, %%rdx\n\t"
      
      "subq %%rax, %2\n\t"
      "subq %%rdx, %3\n\t"
      
      "testq %1, %1\n\t"  // for the branch at the bottom
      
      "xchgq %2, %4\n\t"
      "xchgq %3, %5\n\t"
      
      "jnz 0b\n\t"
      
      : "=r"(t_a), "=r"(t_b), "=r"(u), "=r"(v), "=r"(m), "=r"(n)
      : "0"(t_a), "1"(t_b), "2"(u), "3"(v), "4"(m), "5"(n)
      : "cc", "rax", "rdx");
#else
  int64_t q, t;
  while (b != 0) {
    q = t_a / t_b;
    
    t = b;
    t_b = t_a - q*t_b;
    t_a = t;
    
    t = m;
    m = u - q*m;
    u = t;
    
    t = n;
    n = v - q*n;
    v = t;
  }
#endif

  x = ANTL::negate_using_mask(sm, u);
  y = ANTL::negate_using_mask(sn, v);
  g = t_a;

}