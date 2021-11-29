/**
 * @file qo_cube_plain_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_cube_plain class (for odd char base fields).
 */

template <class T> void CubePlain<T>::cube (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A) {
  T a, b, c, b2, Ca, Cb, Cc;
  T SP, S, v1, u2, v2, N, K, L, TT, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();
  add(b2,b,b);

  // solve SP = v1 (2 b) + u1 a (only need v1)
  XGCD_LEFT (SP, v1, b2, a);

  if (IsOne(SP)) {
    // N = a
    N = a;

    // L = a^2
    sqr(L,a);

    // K = c v1 (v1 (2 b - a c v1) - 2)  (mod L)
    mul(temp,v1,c);
    rem(temp,temp,L);
    mul(K,temp,a);
    rem(K,K,L);
    sub(K,b2,K);
    mul(K,K,v1);
    sub(K,K,2);
    rem(K,K,L);
    mul(K,K,temp);
    rem(K,K,L);
  }
  else {
    // S = u2 (a SP) + v2 (3 b^2 + D)
    mul(SP,SP,a);

    sqr(temp,b);
    mul(temp,temp,3);
    add(temp,temp,Delta);

    XGCD(S,u2,v2,SP,temp);

    // N = a/S
    div(N,a,S);

    // L = N a
    mul(L,N,a);

    // K = -c(v1 u2 a + 2 v2 b)  (mod L)
    mul(K,v1,u2);
    rem(K,K,L);
    mul(K,K,a);
    rem(K,K,L);
    mul(temp,v2,b2);
    rem(temp,temp,L);
    add(K,K,temp);
    mul(K,K,c);
    rem(K,K,L);

    // C = Sc
    mul(c,c,S);
  }

  // T = N K
  mul(TT,N,K);

  // C.a = N L
  mul(Ca,N,L);

  // C.b = b + T
  add(Cb,TT,b);

  // C.c = (S c + K (2 b + T)) / L
  add(Cc,TT,b2);
  //  mul(Cc,Cc,K);
  ::MulExact(Cc,Cc,K,deg(L));
  add(Cc,Cc,c);
  div(Cc,Cc,L);

  C.assign(Ca,Cb,Cc);
}
