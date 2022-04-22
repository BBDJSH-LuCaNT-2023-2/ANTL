/**
 * @file qo_multiply_plain_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_multiply_plain class 
 * (for odd char base fields).
 */

template < class T >  void MultiplyComp<T>::multiply (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B) {
  T a1, a2, b1, b2, c2, Ca, Cb, Cc;
  T SP, S, ab2, v1, u2, v2, K, TT, temp;

  a1 = A.get_a();
  a2 = B.get_a();
  b1 = A.get_b();
  b2 = B.get_b();
  c2 = B.get_c();

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT (SP, v1, a2, a1);

  // K = v1 (b1 - b2) (mod L)
  sub(K,b1,b2);
  mul(K,K,v1);
  rem(K,K,a1);

  if (deg(SP)) {
    add(ab2,b1,b2);

    XGCD(S, u2, v2, SP, ab2);

    // K = u2 K - v2 c2  (mod L)
    mul(K,K,u2);
    mul(temp,v2,c2);
    sub(K,K,temp);

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

  // C.c = (S c2 + K (2 b2 + T)) / L;
  add(Cc,b2,TT);
  add(Cc,Cc,b2);
  //  mul(Cc,Cc,K);
  ::MulExact(Cc,Cc,K,deg(a1));
  add(Cc,Cc,c2);
  div(Cc,Cc,a1);

  C.assign(Ca,Cb,Cc);
  C.reduce();

  //Logic for RelativeGenerator
}
