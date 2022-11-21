/**
 * @file xgcd.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of XGCD_LEFT, XGCD_PARTIAL, and
 * XGCD_PARTIAL_REDUCE
 */

#include <ANTL/XGCD/xgcd.hpp>
#include <ANTL/XGCD/hxgcd.hpp>

//
// XGCD_LEFT functions
//

template <>
void XGCD_LEFT(ZZ & G, ZZ & X, const ZZ & A, const ZZ & B)
{
  ZZ Y;
  NTL::XGCD(G,X,Y,A,B);
}

template <>
void XGCD_LEFT(long & G, long & X, const long & A, const long & B)
{
  long Y;
  NTL::XGCD(G,X,Y,A,B);
}


//
// Partial Euclidean algorithm (for NUCOMP)
//


// Partial Euchlidean algorithm
// Lehmer's version for comuting GCD (for Book's version of NUCOMP, NUDUPL,
//   and NUCUBE algorithm )
//
/*
   Input:  R2 = R_{-1} , R1 = R_{0}, bound
    - R_i is the R - sequence from "Solving the Pell Equation"
      ( R_i = R_{i-2}-q_i R_{i-1} )
    - R1 is reduced modulo R2
   Output: R2 = R_{i-1}, R1 = R_i, C2 = C_{i-1}, C1 = C_i,
    - R_i = 0 or R_i <= bound < R_{i-1}
    - C_i sequence from "Solving the Pell Equation" defined as
      C_{-1}=0, C_{1}=-1  C_i=C_{i-2}-q_i C_{i-1}
   Loop invariant:  abs(C2*R1 - C1*R2) = R_{-1}
*/

void XGCD_PARTIAL(ZZ & R2, ZZ & R1, ZZ & C2, ZZ & C1, const ZZ & bound) {
    ZZ  q, r;

    C2 = 0;
    C1 = -1;

    while (R1 > bound) {
      DivRem(q,r,R2,R1);
      R2 = R1;
      R1 = r;

      r = C2 - q * C1;
      C2 = C1;
      C1 = r;
    }

//   static ZZ q, r, t1, t2;
//   static long A2, A1, TA, B2, B1, TB, rr2, rr1, Tr, qq, bb, T, T1;
//   static int i;
//
//   clear(C1);
//   C2 = to_ZZ(-1);

/*
  ZZ ORIG_R2 = R2;
  ZZ tval = C2*R1 - C1*R2;
  if (abs(tval) != ORIG_R2) {
    cout << "ERROR!" << endl;
    exit(1);
  }
*/

/*
  while (!IsZero(R1) && R1 > bound) {
    T = NumBits (R2) - 31;
    T1 = NumBits (R1) - 31;
    if (T < T1) T=T1;
    if (T < 0) T=0;
    RightShift(r, R2, T);      conv (rr2, r);
    RightShift(r, R1, T);      conv (rr1, r);
    RightShift(r, bound, T);   conv(bb, r);

    A2 = 0;  A1 = 1;
    B2 = 1;  B1 = 0;
    i=0;

    // Euclidean Steps (single precision)
    while ( rr1 != 0  && rr1 > bb ) {
      qq = rr2 / rr1;

      Tr = rr2 - qq*rr1;
      TA = A2 - qq*A1;
      TB = B2 - qq*B1;

      if ( i&1 ) {
        if ( (Tr < -TB) || ( rr1 - Tr < TA - A1 ) ) break;
      }
      else {
        if ( (Tr < -TA) || ( rr1 - Tr < TB - B1 ) ) break;
      }

      rr2 = rr1; rr1 = Tr;
      A2 = A1; A1 = TA;
      B2 = B1; B1 = TB;

      i++;
    }

    if (i == 0) {
      // multiprecsion step
      DivRem(q,R2,R2,R1);
      swap(R2,R1);

      MulSubFrom(C2,C1,q); swap(C2,C1);
    }
    else {
      // recombination
      // r = u*B2 + v*A2;  v = u*B1 + v*A1; u = r;

      mul(r, R2, B2); MulAddTo(r, R1, A2);
      mul(R1, R1, A1); MulAddTo(R1, R2, B1);
      R2 = r;

      // r = p2*A2 + p1*B2;  p2 = p2*A1 + p1*B1; p1 = r;
      mul(r, C2, B2); MulAddTo(r, C1, A2);
      mul(C1, C1, A1); MulAddTo(C1, C2, B1);
      C2 = r;

      if (R1 < 0) { NTL::negate(C1, C1); NTL::negate(R1, R1); }
      if (R2 < 0) { NTL::negate(C2, C2); NTL::negate(R2, R2); }
    }*/

/*
    tval = C2*R1 - C1*R2;
    if (abs(tval) != ORIG_R2) {
      cout << "ERROR!" << endl;
      exit(1);
    }
*/

//   }
//
//   if (R2 < 0) { NTL::negate(C2, C2); NTL::negate(C1, C1); NTL::negate(R2, R2);  }
}

void XGCD_PARTIAL(long & R2, long & R1, long & C2, long & C1, const ZZ & bound) {
  static long q, r, t1, t2;
  static long A2, A1, TA, B2, B1, TB, rr2, rr1, Tr, qq, bb, T, T1;
  static int i;

  clear(C2);
  C1 = -1;

/*
  ZZ ORIG_R2 = R2;
  ZZ tval = C2*R1 - C1*R2;
  if (abs(tval) != ORIG_R2) {
    cout << "ERROR!" << endl;
    exit(1);
  }
*/


  while (!IsZero(R1) && R1 > bound ) {
    T = NumBits (R2) - 31;
    T1 = NumBits (R1) - 31;
    if (T < T1) T=T1;
    if (T < 0) T=0;

    r = R2 >> T;
    //RightShift(r, R2, T);
    conv (rr2, r);

    r = R1 >> T;
    //RightShift(r, R1, T);
    conv (rr1, r);

    r = to_long(bound) >> T;
    //RightShift(r, bound, T);
    //...Not sure if the above is a faithful conversion to long arithmetic, since bound should probably stay as ZZ...
    conv(bb, r);

    A2 = 0;  A1 = 1;
    B2 = 1;  B1 = 0;
    i=0;

    // Euclidean Steps (single precision)
    while ( rr1 != 0  && rr1 > bb ) {
      qq = rr2 / rr1;

      Tr = rr2 - qq*rr1;
      TA = A2 - qq*A1;
      TB = B2 - qq*B1;

      if ( i&1 ) {
        if ( (Tr < -TB) || ( rr1 - Tr < TA - A1 ) ) break;
      }
      else {
        if ( (Tr < -TA) || ( rr1 - Tr < TB - B1 ) ) break;
      }

      rr2 = rr1; rr1 = Tr;
      A2 = A1; A1 = TA;
      B2 = B1; B1 = TB;

      i++;
    }

    if (i == 0) {
      // multiprecsion step

      q = R2/R1;
      R2 = R2 % R1;
      //DivRem(q,R2,R2,R1);
      swap(R2,R1);

      C2 -= C1*q;
      //MulSubFrom(C2,C1,q);
      swap(C2,C1);
    }
    else {
      // recombination
      // r = u*B2 + v*A2;  v = u*B1 + v*A1; u = r;

      mul(r, R2, B2);
      r += R1*A2;
      //MulAddTo(r, R1, A2);

      mul(R1, R1, A1);
      R1 += R2*B1;
      //MulAddTo(R1, R2, B1);
      R2 = r;

      // r = p2*A2 + p1*B2;  p2 = p2*A1 + p1*B1; p1 = r;
      mul(r, C2, B2);
      r += C1*A2;
      //MulAddTo(r, C1, A2);

      mul(C1, C1, A1);
      C1 += C2*B1;
      //MulAddTo(C1, C2, B1);
      C2 = r;

      if (R1 < 0) {
        C1 -= C1;
        R1 -= R1;
//      NTL::negate(C1, C1);
//      NTL::negate(R1, R1);
      }
      if (R2 < 0) {
        C2 -= C2;
        R2 -= C2;
//      NTL::negate(C2, C2);
//      NTL::negate(R2, R2);
      }
    }

/*
    tval = C2*R1 - C1*R2;
    if (abs(tval) != ORIG_R2) {
      cout << "ERROR!" << endl;
      exit(1);
    }
*/

  }

  if (R2 < 0) {
    C2 -= C2;
    C1 -= C2;
    R2 -= R2;
//     NTL::negate(C2, C2);
//     NTL::negate(C1, C1);
//     NTL::negate(R2, R2);
  }
}


//
// Partial Euclidean algorithm (for NUCOMP and fast reduce)
//

void XGCD_PARTIAL(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound)
{
#ifdef TRACE_XGCD
  cout << "XGCD_PARTIAL:" << endl;
#endif
  if (deg(R1) < Thresholds<GF2EX>::get_half_xgcd_partial_crossover(deg(R1),bound))
    XGCD_PARTIAL_ITER(R2,R1,C2,C1,bound);
  else {
    TMatrix<GF2EX> M_out;
    long d_red = deg(R2) - bound;
    HALF_GCD(M_out, R2, R1, d_red);
    clear(C2);
    set(C1);
    mul(C2, C1, M_out);
  }
}



//
// Partial Euclidean algorithm (for fast reduce) with revised termination for
// polynomial types
//

void XGCD_PARTIAL_REDUCE(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even)
{
#ifdef TRACE_XGCD
  cout << "XGCD_PARTIAL_REDUCE:" << endl;
#endif
  if (deg(R1) < Thresholds<GF2EX>::get_half_xgcd_partial_crossover(deg(R1),bound))
    XGCD_PARTIAL_REDUCE_ITER(R2,R1,B2,B1,bound,even);
  else
    HXGCD_PARTIAL_REDUCE(R2,R1,B2,B1,bound,even);
}
