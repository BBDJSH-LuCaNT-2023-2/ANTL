#include <ANTL/Quadratic/Multiply/MultiplyComp.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

template <>
void MultiplyComp<long>::multiply(QuadraticIdealBase<long> &C,
                                const QuadraticIdealBase<long> &A,
                                const QuadraticIdealBase<long> &B) {
  long a1 = A.get_a(), a2 = B.get_a(), b1 = A.get_b(), b2 = B.get_b(),
     c2 = B.get_c();
  long SP, S, ab2, v1, u2, v2, K, T, temp, Ca, Cb, Cc;

  set(*RelativeGenerator);

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT(SP, v1, a2, a1);

  // K = v1 (b1 - b2) / 2 (mod L)
//   sub(K, b1, b2);
//   RightShift(K, K, 1);
//   mul(K, K, v1);
//   rem(K, K, a1);
  K = (v1 * ((b1 - b2) / 2)) % a1;

  if (!IsOne(SP)) {
//     add(ab2, b1, b2);
//     RightShift(ab2, ab2, 1);
    ab2 = (b1 + b2) / 2;

    XGCD(S, u2, v2, SP, ab2);

    // K = u2 K - v2 c2  (mod L)
//     mul(K, K, u2);
//     mul(temp, v2, c2);
//     sub(K, K, temp);
    K = (u2 * K) - (v2 * c2);

    if (!IsOne(S)) {
//       div(a1, a1, S);
      a1 /= S;
//       div(a2, a2, S);
      a2 /= S;
//       mul(c2, c2, S);
      c2 /= S;
    }

//     rem(K, K, a1);
    K %= a1;
  } else {
    S = 1;
  }
  // N = a2;  L = a1;

  // T = NK
//   mul(T, a2, K);
  T = a2 * K;

  // C.a = A.a B.a / d^2 = NL
//   mul(Ca, a2, a1);
  Ca = a2 * a1;

  // C.b = b2 + 2 a2 K = b2 + 2 T
//   LeftShift(Cb, T, 1);
//   add(Cb, Cb, b2);
  Cb = b2 + (2 * T);

  // C.c = (S c2 + K (b2 + T)) / L;
//   add(Cc, b2, T);
//   mul(Cc, Cc, K);
//   add(Cc, Cc, c2);
//   div(Cc, Cc, a1);
  Cc = ((K * (b2 + T)) + c2) / a1;

  C.assign(Ca, Cb, Cc);
  C.reduce();

  ANTL::mul(*RelativeGenerator, *RelativeGenerator,
            *C.get_QO()->get_red_best()->get_RelativeGenerator());

  ANTL::div(*RelativeGenerator, *RelativeGenerator, S);
}
