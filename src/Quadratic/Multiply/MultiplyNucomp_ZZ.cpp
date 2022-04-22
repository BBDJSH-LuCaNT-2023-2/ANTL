/**
 * @file qo_nucomp_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucomp class for ZZ.
 */

#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>
#include <NTL/RR.h>

template <>
void MultiplyNucomp<ZZ>::init(const ZZ &delta_in, const ZZ &h_in, long g_in) {
  MultiplyStrategy<ZZ>::init(delta_in, h_in, 0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}

template <>
void MultiplyNucomp<ZZ>::multiply(QuadraticIdealBase<ZZ> &C,
                                  const QuadraticIdealBase<ZZ> &A,
                                  const QuadraticIdealBase<ZZ> &B) {

  ZZ a1, a2, b1, b2, c2, Ca, Cb, Cc, ss, m;
  ZZ SP, S, v1, u2, v2, K, T, temp;
  ZZ M1, M2, B1, B2, rgA, rgB, rgC, b_temp;

  ZZ TEST, b_mid1, b_mid2, gcd_a1_a2, X, Y, Z, U, R_i, R_i_m1, C_i, C_i_m1,
      floor_root_delta, q, C_temp, R_temp, a_i_p1, b_i_p1, k, B_i, B_i_m1;

  set(*RelativeGenerator);
  // want a1 to be the smaller of the two a coefficients, because initial
  // computations are done mod a1
  if (A.get_a() < B.get_a()) {
    a1 = A.get_a();
    a2 = B.get_a();
    b1 = A.get_b();
    b2 = B.get_b();
    c2 = B.get_c();
  } else {
    a1 = B.get_a();
    a2 = A.get_a();
    b1 = B.get_b();
    b2 = A.get_b();
    c2 = A.get_c();
  }

  // (ss) b_mid1 = (b1 + b2)/2, (m) b_mid2 = (b1 - b2)/2
  b_mid1 = (b1 + b2) / 2;
  b_mid2 = (b1 - b2) / 2;

  // solve (SP) gcd_a1_a2 = X*a2 + v2*a1 (only need X (v1))
  XGCD_LEFT(gcd_a1_a2, X, a2, a1);
  rem(X, X, a1);

  XGCD(S, Y, Z, b_mid1, gcd_a1_a2);

  rem(U, X * Z * b_mid2 - Y * c2, (2 * a1) / S);

  R_i = (2 * a1) / S;
  R_i_m1 = U;

  C_i = 0;
  C_i_m1 = -1;

  Delta = C.get_QO()->get_discriminant();
  floor_root_delta = FloorToZZ(sqrt(to_RR(Delta)));
  NC_BOUND = FloorToZZ(sqrt(to_RR(a1) / to_RR((2 * a2)))* sqrt(sqrt(to_RR(Delta))));

  // check if NUCOMP steps are required
  if (R_i < NC_BOUND) {
    std::cout << "NUCOMP not required!" << std::endl;
    // NUCOMP formulas not required
    a_i_p1 = (a1 * a2) / (S * S);

    b_i_p1 = b2 + (2 * a2 * U) / S;
    rem(b_i_p1, b_i_p1, 2 * a_i_p1);
  }

  else {
    std::cout << "NUCOMP required!" << std::endl;
    // NUCOMP formulas required

    // Execute partial reduction
    //XGCD_PARTIAL(R_i_m1, R_i, C_i_m1, C_i, NC_BOUND);

    // Temporarily using a hard-coded while loop instead of the ought-to-used
    // XGCD_PARTIAL, since it is not currently behaving correctly

    while (R_i >= NC_BOUND) {
      q = R_i_m1 / R_i;

      C_temp = C_i;
      C_i = C_i_m1 - q * C_i;
      C_i_m1 = C_temp;

      R_temp = R_i;
      R_i = R_i_m1 - q * R_i;
      R_i_m1 = R_temp;
    }

//     std::cout << "R_i_m1, R_i, C_i_m1, C_i: " << R_i_m1 << " " << R_i << " "
//               << C_i_m1 << " " << C_i << std::endl;

    M1 = (R_i * a2 + S * (b1 - b2) * C_i) / (2 * a1);
    M2 = (R_i * (b1 + b2) - 2 * S * C_i * c2) / (2*a1 / S);

    if (C_i < 0) {
      a_i_p1 = (R_i * M1 - C_i * M2) / 2;
    }

    else {
      a_i_p1 = (C_i * M2 - R_i * M1) / 2;
    }
    //std::cout << "a_i_p1 is: " << a_i_p1 << std::endl;

    b_i_p1 = (((a2 * R_i + 2 * a_i_p1 * S * C_i_m1) / (S * C_i)) - b2);

    //std::cout << "b_i_p1 is: " << b_i_p1 << std::endl;
  }

  Ca = abs(a_i_p1);

  div(k, floor_root_delta - b_i_p1, 2 * Ca);
  Cb = 2 * Ca * k + b_i_p1;

  Cc = ((Cb * Cb - Delta) / (4 * Ca));

  B_i_m1 = abs(C_i_m1);
  B_i = abs(C_i);

  rgA = S * (2 * a_i_p1 * B_i_m1 + b_i_p1 * B_i);
  rgB = -S * B_i;
  rgC = 2 * a_i_p1;

  RelativeGenerator->set_abd(rgA, rgB, rgC);
  RelativeGenerator->invert();
  if (RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, ZZ(-1));
  }

  // normalize and reduce
  C.assign(Ca, Cb, Cc);
  std::cout << "NUCOMP: Before final reduction, a, b, c are: "
            << Ca << " "
            << Cb << " "
            << Cc << std::endl;
  C.reduce();

  ANTL::mul(*RelativeGenerator, *RelativeGenerator,
            *C.get_QO()->get_red_best()->get_RelativeGenerator());
}
