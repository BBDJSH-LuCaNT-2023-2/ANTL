/**
 * @file qo_nucomp_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucomp class for long.
 */

#include <Quadratic/Multiply/qo_nucomp.hpp>
#include <NTL/RR.h>


template <>
void 
qo_nucomp<long>::init(const long & Din, const long & hin, long gin)
{
  qo_multiply<long>::init(Din,hin,0);
  NC_BOUND = FloorTolong(sqrt(sqrt(abs(to_RR(Delta)))));
}



template <>
void
qo_nucomp<long>::multiply(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A, const QuadraticIdealBase<long> & B)
{
  static long a1, a2, b1, b2, c2, Ca, Cb, Cc, ss, m;
  static long SP, S, v1, u2, v2, K, T, temp;
  static long R1, R2, C1, C2, M1, M2;

  // want a1 to be the smaller of the two a coefficients, because initial
  // computations are done mod a1
  if (A.get_a() < B.get_a()) {
    a1 = A.get_a();
    a2 = B.get_a();
    b1 = A.get_b();
    b2 = B.get_b();
    c2 = B.get_c();
  }
  else {
    a1 = B.get_a();
    a2 = A.get_a();
    b1 = B.get_b();
    b2 = A.get_b();
    c2 = A.get_c();
  }

  // s = (b1 + b2)/2, m = (b1 - b2)/2
  ss = (b1 + b2) >> 1;
  m = (b1 - b2) >> 1;

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT (SP, v1, a2, a1);

  // K = v1 (b1 - b2) / 2 (mod L)
  K = (m*v1) % a1;

  if (SP != 1)
    {
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

  // check if NUCOMP steps are required
  if (a1 < NC_BOUND) {
    // compute with regular multiplication formula (result will be reduced)

  // T = NK
  T = a2*K;

  // C.a = A.a B.a / d^2 = NL
  Ca = a1*a2;

  // C.b = b2 + 2 a2 K = b2 + 2 T
  Cb = (T << 1) + b2;

  // C.c = (S c2 + K (b2 + T)) / L;
  Cc = ((b2+T)*K + c2) / a1;
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // M1 = (N R1 + (b1 - b2) C1 / 2) / L  (T = N R1)
    T = a2*R1;
    M1 = (m*C1 + T) / a1;

    // M2 = (R1(b1 + b2)/2 - c2 S C1) / L
    M2 = (ss*R1 - c2*C1) / a1;

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    if (C1 < 0)
      Ca = R1*M1 - C1*M2;
    else
      Ca = C1*M2 - R1*M1;

    // C.b = 2 (N R1 - C.a C2) / C1 - b2 (mod 2a)
    C.b = ((((T - Ca*C2) << 1) / C1) - b2) % (Ca << 1);

    // C.c = (C.b^2 - Delta) / 4 C.a
    Cc = ((Cb*Cb - Delta) / Ca) >> 2;

    if (Ca < 0) {
      Ca = -Ca;
      Cc = -Cc;
    }
  }

  // normalize
  C.assign(Ca,Cb,Cc);
}
