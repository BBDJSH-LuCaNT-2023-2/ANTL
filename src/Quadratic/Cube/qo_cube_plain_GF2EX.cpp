/**
 * @file qo_cube_plain_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_cube_plain class for GF2EX.
 */

#include <Quadratic/Cube/qo_cube_plain.hpp>

template <> 
void
qo_cube_plain<GF2EX>::cube(QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A)
{
  GF2EX a, b, c, Ca, Cb, Cc;
  GF2EX SP, S, v1, u2, v2, N, K, L, TT, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // solve SP = v1 hx + u1 a (only need v1)
  XGCD_LEFT (SP, v1, hx, a);

  if (IsOne(SP)) {
    // N = a
    N = a;

    // L = a^2
    sqr(L,a);

    // K = c v1^2 (h + a c v1) mod L
    mul(temp,v1,c);
    rem(temp,temp,L);
    mul(K,temp,a);
    rem(K,K,L);
    add(K,K,hx);
    mul(K,K,temp);
    rem(K,K,L);
    mul(K,K,v1);
    rem(K,K,L);
  }
  else {
    // S = u2 (a SP) + v2 (h^2 + QR)
    mul(SP,SP,a);

    sqr(TT,hx);
    mul(temp,a,c);
    add(temp,temp,TT);

    XGCD(S,u2,v2,SP,temp);

    // N = a/S
    div(N,a,S);

    // L = N a = a^2 / S
    mul(L,N,a);

    // K = c(u2 v1 a + v2 h) mod L
    mul(K,u2,v1);
    rem(K,K,L);
    mul(K,K,a);
    rem(K,K,L);
    mul(temp,v2,hx);
    rem(temp,temp,L);
    add(K,K,temp);
    mul(K,K,c);
    rem(K,K,L);

    // c = S c
    mul(c,c,S);
  }

  // T = NK
  mul(TT,N,K);

  // C.a = NL
  mul(Ca,N,L);

  // C.b = b + T
  add(Cb,b,TT);

  //C.c = (S c + K(h + T)) / L
  add(Cc,hx,TT);
  //  mul(Cc,Cc,K);
  ::MulExact(Cc,Cc,K,deg(L));
  add(Cc,Cc,c);
  div(Cc,Cc,L);

  C.assign(Ca,Cb,Cc);
}
