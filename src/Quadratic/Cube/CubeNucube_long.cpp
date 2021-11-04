/**
 * @file qo_nucube_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucube class for long.
  */

#include <ANTL/Quadratic/Cube/CubeNucube.hpp>
#include <NTL/RR.h>


template <> void CubeNucube<long>::init(const long & delta_in, const long & h_in, long g_in) {
  CubeStrategy<long>::init(delta_in,h_in,0);
  SQRT_DELTA = FloorToZZ(sqrt(abs(to_RR(delta_in))));
}


template <> void CubeNucube<long>::cube(QuadraticIdealBase<long> &C, const QuadraticIdealBase<long> &A) {
  static long a, b, c, Ca, Cb, Cc;
  static long SP, S, v1, u2, v2, N, K, L, T, temp, temp2;
  static long B, R1, R2, C1, C2, BB, M1, M2;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // solve SP = v1 b + u1 a (only need v1)
  XGCD_LEFT (SP, v1, b, a);

  if (SP == 1) {
    // N = a
    N = a;

    // L = a^2
    L = a*a;

    // K = c v1 (v1(b - a c v1) - 2) mod L
    v1 %= L;
    temp = (v1*c) % L;

    K = (temp*a) % L;
    K = ((b - K)*v1) % L;
    K = ((K - 2)*temp) % L;
  }
  else {
    // S = u2 (a SP) + v2 (b^2 - ac)
    SP *= a;

    temp = b*b;
    T = a*c;
    temp = temp - T;

    XGCD(S,u2,v2,SP,temp);

    // N = a/S
    N = a/S;

    // L = N a
    L = N*a;

    // K = -c(v1 u2 a + v2 b) mod L
    v2 %= L;
    K = (v1*u2) % L;
    K = (K*a) % L;
    temp = (v2*b) % L;
    K = K + temp;
    K = (K*c) % L;
    K = L - K;

    // C = Sc
    c *= S;
  }

  // Compute NUCOMP termination bound
  B = to_long(SQRT_DELTA) * a;
  B = B >> 1;
  B = SqrRoot(B);
 
  if (L < B) {
    // compute with regular cubing formula (result will be reduced)

    // T = NK
    T = N*K;

    // C.a = N L
    Ca = N*L;

    // C.b = b + 2 T
    Cb = (T << 1) + b;

    // C.c = (S c + K (T + b)) / L
    Cc = ((b+T)*K + c) / L;
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=L; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, ZZ(B));

    // T = N K
    T = (N*K) % L;

    // M1 = (N R1 + T C1) / L  (temp = N R1)
    temp = N*R1;
    M1 = (T*C1 + temp) / L;

    // M2 = (R1(b + T) - c S C1) / L
    M2 = ((b+T)*R1 - c*C1) / L;

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    Ca = R1*M1;
    temp2 = C1*M2;
    if (C1 < 0)
      Ca = Ca - temp2;
    else
      Ca = temp2 - Ca;

    // C.b = 2 (N R1 - C.a C2) / C1 - b (mod 2a)
    Cb = ((((temp - Ca*C2) << 1) / C1) - b) % (Ca << 1);

    // C.c = (C.b^2 - Delta) / 4 C.a
    Cc = ((Cb*Cb - Delta) / Ca) >> 2;

    if (Ca < 0) {
      Ca = -Ca;
      Cc = -Cc;
    }
  }

  C.assign(Ca,Cb,Cc);
}
