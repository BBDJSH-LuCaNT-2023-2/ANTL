/**
 * @file qo_nucube_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucube class for ZZ.
  */

#include <ANTL/Quadratic/Cube/CubeNucube.hpp>
#include <NTL/RR.h>


template <> void CubeNucube<ZZ>::init(const ZZ & delta_in, const ZZ & h_in, long g_in) {
  CubeStrategy<ZZ>::init(delta_in,h_in,0);
  sqrt_delta = FloorToZZ(sqrt(abs(to_RR(delta_in))));
}


template <> void CubeNucube<ZZ>::cube(QuadraticIdealBase<ZZ> &C, const QuadraticIdealBase<ZZ> &A) {
  static ZZ a, b, c, Ca, Cb, Cc;
  static ZZ SP, S, v1, u2, v2, N, K, L, T, temp, temp2;
  static ZZ B, R1, R2, C1, C2, BB, M1, M2;

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

    // K = c v1 (v1(b - a c v1) - 2) mod L
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

    // K = -c(v1 u2 a + v2 b) mod L
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

  // Compute NUCOMP termination bound
  mul(B, sqrt_delta, a);
  RightShift(B,B,1);
  SqrRoot(B, B);
 
  if (L < B) {
    // compute with regular cubing formula (result will be reduced)

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
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=L; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, B);

    // T = N K
    MulMod(T,N,K,L);

    // M1 = (N R1 + T C1) / L  (temp = N R1)
    mul(temp,N,R1);
    mul(M1,T,C1);
    add(M1,M1,temp);
    div(M1,M1,L);

    // M2 = (R1(b + T) - c S C1) / L
    add(M2,b,T);
    mul(M2,M2,R1);
    mul(temp2,c,C1);
    sub(M2,M2,temp2);
    div(M2,M2,L);

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mul(Ca,R1,M1);
    mul(temp2,C1,M2);
    if (sign(C1) < 0)
      sub(Ca,Ca,temp2);
    else
      sub(Ca,temp2,Ca);

    // C.b = 2 (N R1 - C.a C2) / C1 - b (mod 2a)
    mul(Cb,Ca,C2);
    sub(Cb,temp,Cb);
    LeftShift(Cb,Cb,1);
    div(Cb,Cb,C1);
    sub(Cb,Cb,b);
    rem(Cb,Cb,Ca << 1);

    // C.c = (C.b^2 - Delta) / 4 C.a
    sqr(Cc,Cb);
    sub(Cc,Cc,Delta);
    div(Cc,Cc,Ca);
    RightShift(Cc,Cc,2);

    if (Ca < 0) {
      NTL::negate(Ca,Ca);
      NTL::negate(Cc,Cc);
    }
  }

  C.assign(Ca,Cb,Cc);
}
