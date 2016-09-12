/**
 * @file qo_cube_plain_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_cube_plain class for long.
 */

#include <Quadratic/Cube/qo_cube_plain.hpp>

template <> 
void
qo_cube_plain<long>::cube (QuadraticIdealBase<long> &C, const QuadraticIdealBase<long> &A)
{
  static long a, b, c, Ca, Cb, Cc;
  static long SP, S, v1, u2, v2, N, K, L, T, temp;

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
    // T = NK
    T = N*K;

    // C.a = N L
    Ca = N*L;

    // C.b = b + 2 T
    Cb = (T << 1) + b;

    // C.c = (S c + K (T + b)) / L
    Cc = ((T+b)*K + c) / L;

  C.assign(Ca,Cb,Cc);
}
