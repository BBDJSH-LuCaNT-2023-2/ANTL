/**
 * @file xgcd_plain.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of XGCD_LEFT_PLAIN, XGCD_PARTIAL_PLAIN, and
 * XGCD_PARTIAL_REDUCE_PLAIN
 */

#include "../../include/ANTL/XGCD/xgcd_plain.hpp"

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
