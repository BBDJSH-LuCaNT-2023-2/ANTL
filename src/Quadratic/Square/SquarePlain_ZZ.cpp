/**
 * @file qo_square_plain_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_square_plain class for ZZ.
 */

#include <ANTL/Quadratic/Square/SquarePlain.hpp>

template <> void SquarePlain<ZZ>::square (QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A) {
  static ZZ a1, b1, c1, Ca, Cb, Cc;
  static ZZ S, v1, K, T;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();

  // solve S = v1 b1 + u1 a1 (only need v1)
  XGCD_LEFT (S, v1, b1, a1);

  // K = -v1 c1 (mod L)
  mul(K,v1,c1);
  NTL::negate(K,K);

  if (!IsOne (S))
    {
      div(a1,a1,S);
      mul(c1,c1,S);
    }

  rem(K,K,a1);

  // N = L = a1;

  // T = NK
  mul(T,a1,K);

  // C.a = A.a^2 / S^2 = N^2
  sqr(Ca,a1);

  // C.b = b1 + 2 a1 K = b1 + 2 T
  LeftShift(Cb,T,1);
  add(Cb,Cb,b1);

  // C.c = (S c1 + K (b1 + T)) / L;
  add(Cc,b1,T);
  mul(Cc,Cc,K);
  add(Cc,Cc,c1);
  div(Cc,Cc,a1);

  C.assign(Ca,Cb,Cc);
}
