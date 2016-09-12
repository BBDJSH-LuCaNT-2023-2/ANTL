/**
 * @file qo_multiply_plain_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_multiply_plain class for long.
 */

#include <Quadratic/Multiply/qo_multiply_plain.hpp>

template <>
void
qo_multiply_plain<long>::multiply (QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A, const QuadraticIdealBase<long> & B)
{
  static long a1, a2, b1, b2, c2, Ca, Cb, Cc;
  static long SP, S, ab2, v1, u2, v2, K, T, temp;

  a1 = A.get_a();
  a2 = B.get_a();
  b1 = A.get_b();
  b2 = B.get_b();
  c2 = B.get_c();

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT (SP, v1, a2, a1);

  // K = v1 (b1 - b2) / 2 (mod L)
  K = (v1*((b1 - b2) >> 1)) % a1;

  if (SP != 1)
    {
      ab2 = (b1 + b2) >> 1;

      XGCD (S, u2, v2, SP, ab2);

      // K = u2 K - v2 c2  (mod L)
      K = (K*u2 - v2*c2) % a1;

      if (S != 1)
	{
	  a1 /= S;
	  a2 /= S;
	  c2 *= S;
	}
    }

  // N = a2;  L = a1;

  // T = NK
  T = a2*K;

  // C.a = A.a B.a / d^2 = NL
  Ca = a1*a2;

  // C.b = b2 + 2 a2 K = b2 + 2 T
  Cb = (T << 1) + b2;

  // C.c = (S c2 + K (b2 + T)) / L;
  Cc = ((b2+T)*K + c2) / a1;

  C.assign(Ca,Cb,Cc);
}

