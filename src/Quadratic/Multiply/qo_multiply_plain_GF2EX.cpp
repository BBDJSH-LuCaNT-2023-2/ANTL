/**
 * @file qo_multiply_plain_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_multiply_plain class for GF2EX.
 */

#include <Quadratic/Multiply/qo_multiply_plain.hpp>

template <> 
void
qo_multiply_plain<GF2EX>::multiply (QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A, const QuadraticIdealBase<GF2EX> &B)
{
  GF2EX a1, a2, b1, b2, c2, Ca, Cb, Cc;
  GF2EX SP, S, ab2, v1, u2, v2, K, TT, temp;

  a1 = A.get_a();
  a2 = B.get_a();
  b1 = A.get_b();
  b2 = B.get_b();
  c2 = B.get_c();

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT (SP, v1, a2, a1);

  // K = v1 (b1 + b2) (mod L)
  add(ab2,b1,b2);
  mul(K,ab2,v1);
  rem(K,K,a1);

  if (deg(SP)) {
    add(ab2,ab2,hx);

    XGCD(S, u2, v2, SP, ab2);

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

  // T = NK
  mul(TT,a2,K);

  // C.a = A.a B.a / d^2 = NL
  mul(Ca,a2,a1);

  // C.b = b2 + a2 K = b2 + T
  add(Cb,TT,b2);

  // C.c = (S c2 + K (hx + T)) / L;
  add(Cc,hx,TT);
  //  mul(Cc,Cc,K);
  ::MulExact(Cc,Cc,K,deg(a1));
  add(Cc,Cc,c2);
  div(Cc,Cc,a1);

  C.assign(Ca,Cb,Cc);
}
