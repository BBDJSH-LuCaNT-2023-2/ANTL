/**
 * @file qo_square_plain_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_square_plain class (for odd char base fields).
 */

template <class T> void SquarePlain<T>::square (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A) {
  T a1, b1, c1, b2, Ca, Cb, Cc;
  T S, v1, K, TT;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();
  add(b2,b1,b1);

  // solve S = v1 (2 b1) + u1 a1 (only need v1)
  XGCD_LEFT (S, v1, b2, a1);

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
  mul(TT,a1,K);

  // C.a = A.a^2 / S^2 = N^2
  sqr(Ca,a1);

  // C.b = b1 + N K = b1 + T
  add(Cb,TT,b1);

  // C.c = (S c1 + K (2 b1 + T)) / L;
  add(Cc,b2,TT);
  //  mul(Cc,Cc,K);
  ::MulExact(Cc,Cc,K,deg(a1));
  add(Cc,Cc,c1);
  div(Cc,Cc,a1);

  C.assign(Ca,Cb,Cc);
}
