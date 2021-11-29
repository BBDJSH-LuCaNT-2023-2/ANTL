/**
 * @file qo_square_plain_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_square_plain class for long.
 */

#include <ANTL/Quadratic/Square/SquarePlain.hpp>

template <> void SquarePlain<long>::square (QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A) {
  static long a1, b1, c1, Ca, Cb, Cc;
  static long S, v1, K, T;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();

  // solve S = v1 b1 + u1 a1 (only need v1)
  XGCD_LEFT (S, v1, b1, a1);

  // K = -v1 c1 (mod L)
  K = (-v1*c1) % a1;

  if (S != 1)
    {
      a1 /= S;
      c1 /= S;
    }

  // N = L = a1;

  // T = NK
  T = a1*K;

  // C.a = A.a^2 / S^2 = N^2
  Ca = a1*a1;

  // C.b = b1 + 2 a1 K = b1 + 2 T
  Cb = (T << 1) + b1;

  // C.c = (S c1 + K (b1 + T)) / L;
  Cc = ((b1+T)*K + c1) / a1;

  C.assign(Ca,Cb,Cc);
}
