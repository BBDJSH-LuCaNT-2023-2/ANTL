/**
 * @file qo_nudupl_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nudupl class for GF2EX.
 */

#include <Quadratic/Square/qo_nudupl.hpp>

template <> 
void
qo_nudupl<GF2EX>::square(QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A)
{
  GF2EX a1, b1, c1, Ca, Cb, Cc;
  GF2EX S, v1, K, TT;
  GF2EX R1, R2, C1, C2, M2, temp;
  long BOUND;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();

  // solve S = v1 hx + u1 a1 (only need v1)
  XGCD_LEFT (S, v1, hx, a1);

  // K = v1 c1 (mod L)
  mul(K,v1,c1);

  if (!IsOne (S))
    {
      div(a1,a1,S);
      mul(c1,c1,S);
    }

  rem(K,K,a1);

  // N = L = a1

  // check if NUCOMP steps are required
  if ((deg(a1) << 1) <= genus) {
    // compute with regular multiplication formula (result will be reduced)

    // T = NK
    mul(TT,a1,K);

    // C.a = A.a^2 / S^2 = N^2
    sqr(Ca,a1);

    // C.b = b1 + N K = b1 + T
    add(Cb,TT,b1);

    // C.c = (S c1 + K (hx + T)) / L;
    add(Cc,hx,TT);
    //    mul(C.c,C.c,K);
    ::MulExact(Cc,Cc,K,deg(a1));
    add(Cc,Cc,c1);
    div(Cc,Cc,a1);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    if (deg(b1) == genus + 1)
      BOUND = (genus + 1) >> 1;
    else
      BOUND = (genus) >> 1;

    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, BOUND);

    // M1 = R1

    // M2 = (R1 hx - c1 S C1) / L
    //   mul(M2,qo<GF2EX>::hx,R1);
    //    mul(temp,c1,C1);
    ::MulExact(M2,hx,R1,deg(a1));
    ::MulExact(temp,c1,C1,deg(a1));
    add(M2,M2,temp);
    div(M2,M2,a1);

    // C.a = R1^2 + C1 M2
    sqr(Ca,R1);
    mul(temp,C1,M2);
    add(Ca,Ca,temp);

    // C.b = (N R1 + C.a C2) / C1 + b1 + hx (mod a)
    //   mul(Cb,Ca,C2);
    //   mul(temp,a1,R1);
    ::MulExact(Cb,Ca,C2,deg(C1));
    ::MulExact(temp,a1,R1,deg(C1));
    add(Cb,temp,Cb);
    div(Cb,Cb,C1);
    add(Cb,Cb,b1);
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
