/**
 * @file qo_cube_plain_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_cube_plain class for ZZ.
 */

#include <Quadratic/Cube/qo_cube_plain.hpp>

template <> 
void
qo_cube_plain<ZZ>::cube (QuadraticIdealBase<ZZ> &C, const QuadraticIdealBase<ZZ> &A)
{
  static ZZ a, b, c, Ca, Cb, Cc;
  static ZZ SP, S, v1, u2, v2, N, K, L, T, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // solve SP = v1 b + u1 a (only need v1)
  XGCD_LEFT (SP, v1, b, a);

  if (IsOne(SP)) {
    // N = a
    N = a;

    // L = a^2
    sqr(L,a);

    // K = c v1 (v1(b - a c v1) - 2) (mod L)
    rem(v1,v1,L);
    mul(temp,v1,c);
    rem(temp,temp,L);

    MulMod(K,temp,a,L);
    SubMod(K,b,K,L);
    MulMod(K,K,v1,L);
    SubMod(K,K,2,L);
    MulMod(K,K,temp,L);
  }
  else {
    // S = u2 (a SP) + v2 (b^2 - ac)
    mul(SP,SP,a);

    sqr(temp,b);
    mul(T,a,c);
    sub(temp,temp,T);

    XGCD(S,u2,v2,SP,temp);

    // N = a/S
    div(N,a,S);

    // L = N a
    mul(L,N,a);

    // K = -c(v1 u2 a + v2 b) (mod L)
    rem(v2,v2,L);
    mul(K,v1,u2);
    rem(K,K,L);
    MulMod(K,K,a,L);
    MulMod(temp,v2,b,L);
    AddMod(K,K,temp,L);
    MulMod(K,K,c,L);
    sub(K,L,K);

    // C = Sc
    mul(c,c,S);
  }

  // T = NK
  mul(T,N,K);

  // C.a = N L
  mul(Ca,N,L);

  // C.b = b + 2 T
  LeftShift(Cb,T,1);
  add(Cb,Cb,b);

  // C.c = (S c + K (T + b)) / L
  add(Cc,T,b);
  mul(Cc,Cc,K);
  add(Cc,Cc,c);
  div(Cc,Cc,L);

  C.assign(Ca,Cb,Cc);
}
