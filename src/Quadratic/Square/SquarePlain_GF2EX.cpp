/**
 * @file qo_square_plain_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_square_plain class for GF2EX.
 */

#include <ANTL/Quadratic/Square/SquarePlain.hpp>

template <> void SquarePlain<GF2EX>::square (QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A) {
  GF2EX a1, b1, c1, Ca, Cb, Cc;
  GF2EX S, v1, K, TT;

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

  // N = L = a1;

  // T = NK
  mul(TT,a1,K);

  // C.a = A.a^2 / S^2 = N^2
  sqr(Ca,a1);

  // C.b = b1 + N K = b1 + T
  add(Cb,TT,b1);

  // C.c = (S c1 + K (hx + T)) / L;
  add(Cc,hx,TT);
  //  mul(Cc,Cc,K);
  ::MulExact(Cc,Cc,K,deg(a1));
  add(Cc,Cc,c1);
  div(Cc,Cc,a1);

  C.assign(Ca,Cb,Cc);
}
