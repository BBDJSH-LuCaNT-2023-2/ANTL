/**
 * @file qo_nudupl_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nudupl class for long.
 */

#include <ANTL/Quadratic/Square/SquareNudupl.hpp>
#include <NTL/RR.h>


template <> void SquareNudupl<long>::init(const long & delta_in, const long & h_in, long g_in) {
  SquareStrategy<long>::init(delta_in,h_in,0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}


template <> void SquareNudupl<long>::square(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A) {
  static long a1, b1, c1, Ca, Cb, Cc;
  static long S, v1, K, T, temp;
  static long R1, R2, C1, C2, M2;

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

  // check if NUCOMP steps are required
  if (a1 < NC_BOUND) {
    // compute with regular squaring formula (result will be reduced)

  // T = NK
  T = a1*K;

  // C.a = A.a^2 / S^2 = N^2
  Ca = a1*a1;

  // C.b = b1 + 2 a1 K = b1 + 2 T
  Cb = (T << 1) + b1;

  // C.c = (S c1 + K (b1 + T)) / L;
  Cc = ((b1+T)*K + S*c1) / a1;
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // M1 = R1

    // M2 = (R1 b1 - c1 S C1) / L
    M2 = (R1*b1 - c1*C1) / a1;

    // C.a = (-1)^(i-1) (R1^2 - C1 M2)
    if (C1 < 0)
      Ca = R1*R1 - C1*M2;
    else
      Ca = C1*M2 - R1*R1;

    // C.b = 2 (N R1 - C.a C2) / C1 - b1 (mod 2a)
    Cb = ((((a1*R1 - Ca*C2) << 1) / C1) - b1) % (Ca << 1);

    // C.c = (C.b^2 - Delta) / 4 C.a
    Cc = ((Cb*Cb - Delta) / Ca) >> 2;

    if (Ca < 0) {
      Ca = -Ca;
      Cc = -Cc;
    }
  }

  // normalize and reduce
  C.assign(Ca,Cb,Cc);
}
