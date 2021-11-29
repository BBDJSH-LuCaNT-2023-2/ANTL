/**
 * @file qo_nucomp_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucomp class for GF2EX.
 */

#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>

template <> void MultiplyNucomp<GF2EX>::multiply(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A, const QuadraticIdealBase<GF2EX> & B) {
  GF2EX a1, a2, b1, b2, c2, Ca, Cb, Cc, ss, m;
  GF2EX SP, S, v1, u2, v2, K, TT, temp;
  GF2EX R1, R2, C1, C2, M1, M2;
  long BOUND;

  // want a1 to have the smaller degree, because initial computations are
  // done mod a1
  if (deg(A.get_a()) < deg(B.get_a())) {
    a1 = A.get_a();
    a2 = B.get_a();
    b1 = A.get_b();
    b2 = B.get_b();
    c2 = B.get_c();
  }
  else {
    a1 = B.get_a();
    a2 = A.get_a();
    b1 = B.get_b();
    b2 = A.get_b();
    c2 = A.get_c();
  }

  // s = b1+b2
  add(ss,b1,b2);
  add(m,ss,hx);

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT(SP, v1, a2, a1);

  // K = v1 (b1 + b2) (mod L)
  mul(K,ss,v1);
  rem(K,K,a1);

  if (deg(SP)) {
    XGCD(S, u2, v2, SP, m);

    // K = u2 K + v2 c2  (mod L)
    mul(K,K,u2);
    mul(temp,v2,c2);
    add(K,K,temp);

    if (!IsOne(S)) {
      div(a1,a1,S);
      div(a2,a2,S);
      mul(c2,c2,S);
    }

    rem(K,K,a1);
  }

  // N = a2;  L = a1;

  // check if NUCOMP steps are required
  if (deg(a1) + deg(a2) <= genus) {
    // compute with regular multiplication formula (result will be reduced)

    // T = NK
    mul(TT,a2,K);

    // C.a = A.a B.a / d^2 = NL
    mul(Ca,a2,a1);

    // C.b = b2 + a2 K = b2 + T
    add(Cb,TT,b2);

    // C.c = (S c2 + K (hx + T)) / L;
    add(Cc,hx,TT);
    //    mul(Cc,Cc,K);
    ::MulExact(Cc,Cc,K,deg(a1));
    add(Cc,Cc,c2);
    div(Cc,Cc,a1);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    if (deg(b2) == genus + 1)
      BOUND = (deg(a1) - deg(a2) + genus + 1) >> 1;
    else
      BOUND = (deg(a1) - deg(a2) + genus) >> 1;

    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, BOUND);

    // M1 = (N R1 + (b1 + b2) C1) / L  (T = N R1)
    mul(TT,a2,R1);
    //    mul(M1,ss,C1);
    ::MulExact(M1,ss,C1,deg(a1));
    add(M1,M1,TT);
    div(M1,M1,a1);

    // M2 = (R1(b1 + b2 + hx) + c2 S C1) / L
    //    mul(M2,m,R1);
    //   mul(temp,c2,C1);
    ::MulExact(M2,m,R1,deg(a1));
    ::MulExact(temp,c2,C1,deg(a1));
    add(M2,M2,temp);
    div(M2,M2,a1);

    // C.a = R1 M1 + C1 M2
    mul(Ca,R1,M1);
    mul(temp,C1,M2);
    add(Ca,Ca,temp);

    // C.b = (N R1 + C.a C2) / C1 + b2 + hx (mod a)
    //   mul(Cb,Ca,C2);
    ::MulExact(Cb,Ca,C2,deg(C1));
    add(Cb,TT,Cb);
    div(Cb,Cb,C1);
    add(Cb,Cb,b2);
    add(Cb,Cb,hx);
    rem(Cb,Cb,Ca);

    // C.c = (C.b^2 + C.b hx + Delta) / C.a
    add(Cc,Cb,hx);
    //    mul(Cc,Cc,Cb);
    ::MulExact(Cc,Cc,Cb,deg(Ca));
    add(Cc,Cc,Delta);
    div(Cc,Cc,Ca);
  }

  // normalize and reduce
  C.assign(Ca,Cb,Cc);
  C.reduce();
}
