/**
 * @file qo_nucube_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucube class for GF2EX.
 */

#include <Quadratic/Cube/qo_nucube.hpp>

template <> 
void
qo_nucube<GF2EX>::cube(QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A)
{
  GF2EX a, b, c, Ca, Cb, Cc;
  GF2EX SP, S, v1, u2, v2, N, K, L, TT;
  GF2EX B, R1, R2, C1, C2, BB, M1, M2, temp, temp2;
  long BOUND;

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

  // Compute NUCOMP termination bound
  BOUND = (deg(a) + genus) >> 1;

  // check if NUCOMP steps are required
  if (deg(L) <= BOUND) {
    // compute with regular multiplication formula (result will be reduced)

    // T = N K
    mul(TT,N,K);

    // C.a = N L
    mul(Ca,N,L);

    // C.b = b + T
    add(Cb,TT,b);

    // C.c = (S c + K (hx + T)) / L
    add(Cc,hx,TT);
    //    mul(Cc,Cc,K);
    ::MulExact(Cc,Cc,K,deg(L));
    add(Cc,Cc,c);
    div(Cc,Cc,L);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=L; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, BOUND);

    // T = N K
    MulMod(TT,N,K,L);

    // M1 = (N R1 + TT C1) / L  (T = N R1)
    mul(temp,N,R1);
    //    mul(M1,TT,C1);
    ::MulExact(M1,TT,C1,deg(L));
    add(M1,M1,temp);
    div(M1,M1,L);

    // M2 = (R1 (hx + TT) + c2 S C1) / L
    add(M2,TT,hx);
    //   mul(M2,M2,R1);
    //    mul(temp2,c,C1);
    ::MulExact(M2,M2,R1,deg(L));
    ::MulExact(temp2,c,C1,deg(L));
    add(M2,M2,temp2);
    div(M2,M2,L);

    // C.a = R1 M1 + C1 M2
    mul(Ca,R1,M1);
    mul(temp2,C1,M2);
    add(Ca,Ca,temp2);

    // C.b = (N R1 + C.a C2) / C1 + b + hx (mod a)
    //    mul(Cb,Ca,C2);
    ::MulExact(Cb,Ca,C2,deg(C1));
    add(Cb,Cb,temp);
    div(Cb,Cb,C1);
    add(Cb,Cb,b);
    add(Cb,Cb,hx);
    rem(Cb,Cb,Ca);

    // C.c = (C.b^2 + C.b hx + Delta) / C.a
    add(Cc,Cb,hx);
    //    mul(C.c,C.c,C.b);
    ::MulExact(Cc,Cc,Cb,deg(Ca));
    add(Cc,Cc,Delta);
    div(Cc,Cc,Ca);
  }

  // normalize and reduce
  C.assign(Ca,Cb,Cc);
  C.reduce();
}

