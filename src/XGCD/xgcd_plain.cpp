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
//Plain XGCD from liboptarith 
//Written by Maxwell Sayles 
template<> 
void XGCD_PLAIN(long & G, long & X, long & Y, const long & A, const long & B){
//   std::cout << "using XGCD_PLAIN" << std::endl;

  long t_a = A;
  long t_b = B;

  ulong sm = A >> 63;
  ulong sn = B >> 63;
  t_a = ANTL::negate_using_mask<long>(sm, A);
  t_b = ANTL::negate_using_mask<long>(sn, B);
    
  long m = 0;
  long n = 1;
  long u = 1;
  long v = 0;
  
  if (t_a == 0) {
    X = 1;
    Y = 0;
    G = B;
    return;
  }
  if (t_b == 0) {
    X = 0;
    Y = 1;
    G = A;
    return;
  }

#if defined(__x86_64)
//     std::cout << "using XGCD_PLAIN ASM" << std::endl;
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
//   std::cout << "using XGCD_PLAIN NON-ASM" << std::endl;
  long q, t;
  while (B != 0) {
    q = t_a / t_b;
    
    t = B;
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

  X = ANTL::negate_using_mask<long>(sm, u);
  Y = ANTL::negate_using_mask<long>(sn, v);
  G = t_a;

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

template<>
void XGCD_LEFT_PLAIN(long & G, long & X, const long & A, const long & B) {
//   std::cout << "using XGCD_LEFT_PLAIN" << std::endl;
  long ma = A;
  long mb = B;

  long ta;
  int sa;
#if !defined(__x86_64)
  long q, t;
#endif

  if (ma < 0) {
    sa = -1;
    ma = -ma;
  } else {
    sa = 1;
  }
  if (mb < 0) {
    mb = -mb;
  }

  ta = 0;
  X = 1;

  if (mb == 0) {
    G = ma;
    return;
  }
  if (ma == 0) {
    X = 0;
    G = mb;
    return;
  }

#if defined(__x86_64)
//   std::cout << "using XGCD_LEFT_PLAIN ASM" << std::endl;
  // 64bit gcd
  asm("0:\n\t"
      "xorq %%rdx, %%rdx\n\t"
      "movq %1, %%rax\n\t"
      "divq %2\n\t"           // rdx = a%b, rax = a/b
      "movq %2, %1\n\t"       // a = b
      "movq %%rdx, %2\n\t"    // b = rdx

      "imulq %3, %%rax\n\t"   // rax = q*a
      "subq %%rax, %0\n\t"    // u -= q*a
      "movq %3, %%rax\n\t"
      "movq %0, %3\n\t"
      "movq %%rax, %0\n\t"

      "testq %2, %2\n\t"
      "jz 1f\n\t"
      "cmpq %4, %1\n\t"       // 64bit constants not permitted in compare
      "jg 0b\n\t"             // unsigned (ja) did not work
      "cmpq %4, %0\n\t"
      "jg 0b\n\t"

      "1:\n\t"
      : "=r"(X), "=r"(ma), "=r"(mb), "=r"(ta)
      : "0"(X), "1"(ma), "2"(mb), "3"(ta), "r"((long)0x7FFFFFFFLL)
      : "cc", "rax", "rdx");
  // either b == 0, or both a and b are 32bit

  if (mb != 0) {
    uint32_t a32 = ma;
    uint32_t b32 = mb;
    // 32bit gcd
    asm("0:\n\t"
	"xorl %%edx, %%edx\n\t"
	"movl %1, %%eax\n\t"
	"divl %2\n\t"           // rdx = a%b, rax = a/b
	"movl %2, %1\n\t"       // a = b
	"movl %%edx, %2\n\t"    // b = rdx

	"imulq %3, %%rax\n\t"   // rax = q*a
	"subq %%rax, %0\n\t"    // u -= q*a
	"movq %3, %%rax\n\t"
	"movq %0, %3\n\t"
	"movq %%rax, %0\n\t"

	"testl %2, %2\n\t"
	"jnz 0b\n\t"
	: "=r"(X), "=r"(a32), "=r"(b32), "=r"(ta)
	: "0"(X), "1"(a32), "2"(b32), "3"(ta)
	: "cc", "rax", "rdx");
    ma = a32;
    mb = b32;
  }
#else
//   std::cout << "using XGCD_LEFT_PLAIN NON-ASM" << std::endl;
  while (B != 0) {
    q = A / B;

    t = B;
    B = A - q*B;
    A = t;

    t = ta;
    ta = (X) - q*ta;
    X = t;
  }
#endif
  X *= sa;
  G = ma;
}

template<>
void XGCD_LEFT_PLAIN(ZZ & G, ZZ & X, const ZZ & A, const ZZ & B){
  //A must be >= B
  mpz_t a,b,g,x;
  assert(A>=B);
  mpz_init_set_ui(a,0);
  mpz_init_set_ui(b,0);

  //ZZ_limbs_get only returns a readonly pointer
  //must convert to mpz_ts and get writable pointers
  ANTL::ZZToMpz(A,a);
  ANTL::ZZToMpz(B,b);
  mpz_init_set(g,a);
  mpz_init_set(x,a);

  mp_size_t a_size, b_size, g_size, x_size;
  a_size = A.size();
  b_size = B.size();
  mp_limb_t *g_limbs, *a_limbs, *b_limbs, *x_limbs;
  a_limbs = mpz_limbs_write(a, a_size);
  b_limbs = mpz_limbs_write(b, b_size);
  g_limbs = mpz_limbs_write(g, a_size);
  x_limbs = mpz_limbs_write(x, a_size);
  g_size = mpn_gcdext(g_limbs, x_limbs, &x_size, a_limbs, a_size, b_limbs, b_size);
  //correct sign on x_size
  cout << "x_size premod: " << x_size << endl;
  int sign = x_size/abs(x_size);
  x_size=abs(x_size);
  cout << "x_size postmod: " << x_size << endl;
  //rebuild ZZ values
  ZZ_limbs_set(G, g_limbs, g_size);
  ZZ_limbs_set(X, x_limbs, x_size);
  X = sign*X;

  //(g-ax)/b
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
