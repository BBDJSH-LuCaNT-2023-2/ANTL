#include <ANTL/Quadratic/Multiply/MultiplyNucomp_Opt.hpp>
#include <NTL/RR.h>

void temp_invert_and_normalize_nucomp(long &rel_gen_a, long &rel_gen_b, long &rel_gen_d, long &delta);

template <>
void MultiplyNucompOpt<long>::construct_relative_generator(
    long &rel_gen_a, long &rel_gen_b, long &rel_gen_d, QuadraticIdealBase<long> &C,
    long OB, long BB, long S);

template <>
void MultiplyNucompOpt<long>::init(const long &delta_in, const long &h_in,
                                 long g_in) {
  init(delta_in, h_in, 0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}

template <>
void MultiplyNucompOpt<long>::multiply(QuadraticIdealBase<long> &C,
                                     const QuadraticIdealBase<long> &A,
                                     const QuadraticIdealBase<long> &B) {

  Delta = C.get_QO()->get_discriminant();
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));

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
  } else {
    a1 = B.get_a();
    a2 = A.get_a();
    b1 = B.get_b();
    b2 = A.get_b();
    c2 = A.get_c();
  }

  // s = (b1 + b2)/2, m = (b1 - b2)/2
  ss = (b1 + b2) / 2;

  m = (b1 - b2) / 2;

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT(SP, v1, a2, a1);

  // K = v1 (b1 - b2) / 2 (mod L)
//   std::cout << "NUCOMP: Before K = (m * v1) mod a1" << std::endl;
//   std::cout << "NUCOMP: m is " << m << std::endl;
//   std::cout << "NUCOMP: v1 is " << v1 << std::endl;
//   std::cout << "NUCOMP: a1 is " << a1 << std::endl;
  K = (m * v1) % a1;
//   std::cout << "NUCOMP: K is " << K << std::endl;

  S = 1;
  if (!IsOne(SP)) {
    XGCD(S, u2, v2, SP, ss);

    // K = u2 K - v2 c2 (mod L)
    K *= u2;
    K -= v2*c2;
//     std::cout << "NUCOMP: K is " << K << std::endl;

    if (!IsOne(S)) {
      a1 /= S;
      a2 /= S;
      c2 *= S;
    }

    K %= a1;
//     std::cout << "NUCOMP: K is " << K << std::endl;
  }

  if (K < 0) {
    K += a1;
  }

  // N = a2;  L = a1;

  // check if NUCOMP steps are required
  if (a1 <= NC_BOUND) {
    // compute with regular multiplication formula (result will be reduced)

    // T = NK
    T = a2 * K;

    // C.a = A.a B.a / d^2 = NL
    Ca = a2 * a1;

    // C.b = b2 + 2 a2 K = b2 + 2 T
    Cb = b2 + (2 * T);

    // C.c = (S c2 + K (b2 + T)) / L;
//     std::cout << "NUCOMP: b2 is " << b2 << std::endl;
//     std::cout << "NUCOMP: T is " << T << std::endl;
//     std::cout << "NUCOMP: K is " << K << std::endl;
//     std::cout << "NUCOMP: c2 is " << c2 << std::endl;
//     std::cout << "NUCOMP: a1 is " << a1 << std::endl;
    Cc = (c2 + (K * (b2 + T))) / a1;

    // Set a, b, c (DO NOT ASSIGN/NORMALIZE)
//     std::cout << "NUCOMP: Ca, Cb, Cc are " << Ca << ", " << Cb << ", " << Cc << std::endl;
    C.set_a(Ca);
    C.set_b(Cb);
    C.set_c(Cc);

    C2 = 1;
    C1 = 0;
  } else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2 = a1;
    R1 = K;

//     std::cout << "NUCOMP: Before XGCD_PARTIAL" << std::endl;
//     std::cout << "NUCOMP: R2 is " << R2 << std::endl;
//     std::cout << "NUCOMP: R1 is " << R1 << std::endl;
//     std::cout << "NUCOMP: C2 is " << C2 << std::endl;
//     std::cout << "NUCOMP: C1 is " << C1 << std::endl;

    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

//     std::cout << "NUCOMP: After XGCD_PARTIAL" << std::endl;
//     std::cout << "NUCOMP: R2 is " << R2 << std::endl;
//     std::cout << "NUCOMP: R1 is " << R1 << std::endl;
//     std::cout << "NUCOMP: C2 is " << C2 << std::endl;
//     std::cout << "NUCOMP: C1 is " << C1 << std::endl;

    // M1 = (N R1 + (b1 - b2) C1 / 2) / L  (T = N R1)
    M1 = ((m * C1) + (a2 * R1)) / a1;

    // M2 = (R1(b1 + b2)/2 - c2 S C1) / L
    M2 = ((ss * R1) - (c2 * C1)) / a1;
//     std::cout << "NUCOMP: M1 is " << M1 << std::endl;
//     std::cout << "NUCOMP: M2 is " << M2 << std::endl;

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)

    if (C1 > 0) {
      Ca = (R1 * M1) - (C1 * M2);
    }
    else {
      Ca = (C1 * M2) - (R1 * M1);
    }

    // C.b = 2 (N R1 - C.a C2) / C1 - b1
    Cb = (a2 * R1 + C2 * Ca) << 1;
    Cb /= C1;
    Cb -= b2;

    // C.c = (C.b^2 - Delta) / 4 C.a
    Cc = (((Cb * Cb) - Delta) / Ca) / 4;


    // Set a, b, c (DO NOT ASSIGN/NORMALIZE)
//     std::cout << "NUCOMP: Ca, Cb, Cc are " << Ca << ", " << Cb << ", " << Cc << std::endl;
    C.set_a(Ca);
    C.set_b(Cb);
    C.set_c(Cc);

    // Partial update of the distance
    // arb_add(qie->distance, qie->distance, qie->distance, precision);
  }

  // Reduce and get the coefficients of the relative generator
//   std::cout << "NUCOMP: Before C.reduce(), C is" << C << std::endl;
  C.reduce();
//   std::cout << "NUCOMP: After C.reduce(), C is" << C << std::endl;

  // relative_generator
  long rel_gen_a, rel_gen_b, rel_gen_d;
  RR relative_generator;

  construct_relative_generator(rel_gen_a, rel_gen_b, rel_gen_d, C, abs(C2),
                               abs(C1), S);

  temp_invert_and_normalize_nucomp(rel_gen_a, rel_gen_b, rel_gen_d, Delta);

  RelativeGenerator->set_abd(rel_gen_a, rel_gen_b, rel_gen_d);
  if (RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, -1);
  }

}

//
// construct relative generator
// computes the coefficients a, b, d of the relative generator
// gamma = (a + b*sqrt(Delta))/d such that A^2 = (gamma) B
// Pass OB=1 and BB=0
//
template <>
void MultiplyNucompOpt<long>::construct_relative_generator(
    long &rel_gen_a, long &rel_gen_b, long &rel_gen_d, QuadraticIdealBase<long> &C,
    long OB, long BB, long S) {
  static long NB;
  bool con_rel_gen_dbg = false;

  if (con_rel_gen_dbg) {
    // printf("\n-->--> construct_relative_generator:\n");
    // printf("OB=%ld, BB=%ld\n",OB,BB);
    std::cout << "\n-->--> construct_relative_generator:" << std::endl;
    std::cout << "OB=" << OB << ", BB=" << BB << std::endl;
  }

  for (long i = 0; i < C.get_num_q(); ++i) {
    NB = C.get_qlist_i(i) * BB + OB;
    OB = BB;
    BB = NB;

    if (con_rel_gen_dbg) {
      // printf("&&& q=%ld, B=%ld\n",qie->qlist[i],BB);
      std::cout << "&&& q=" << C.get_qlist_i(i) << ", B=" << BB << std::endl;
    }
  }

  rel_gen_a = S * ((OB * C.get_a() << 1) + BB * C.get_b());
  rel_gen_b = -S * BB;
  rel_gen_d = (C.get_a() << 1);

  if (con_rel_gen_dbg) {
    // printf("OB=%ld, BB=%ld, a=%ld, b=%ld, S=%ld\n",OB,BB,qie->a,qie->b,S);
    // printf("relgen_a=%ld, relgen_b=%ld,
    // relgen_d=%ld\n",*rel_gen_a,*rel_gen_b,*rel_gen_d); printf("-->--> done
    // construct_relative_generator\n\n");

    std::cout << "OB=" << OB << ", BB=" << BB << ", a=" << C.get_a()
              << ", b=" << C.get_b() << ", S=" << S << std::endl;
    std::cout << "relgen_a=" << rel_gen_a << ", relgen_b=" << rel_gen_b
              << ", relgen_d=" << rel_gen_d << std::endl;
    std::cout << "-->--> done construct_relative_generator\n\n" << std::endl;
  }
}

void temp_invert_and_normalize_nucomp(long &rel_gen_a, long &rel_gen_b, long &rel_gen_d, long &delta) {

  ZZ a, b, d, newA, newB, newD, temp;

  a = rel_gen_a;
  b = rel_gen_b;
  d = rel_gen_d;

  // ((a + b rho) / d)^-1 = (ad - bd rho) / (a^2 - b^2 Delta)
  newA = a * d;

  newB = b * d;

  newD =  a * a;
  temp =  b * b;
  temp = temp * delta;
  newD = newD - temp;

  a = newA;
  b = -newB;
  d = newD;

  if (d < 0) {
    NTL::negate(a, a);
    NTL::negate(b, b);
    NTL::negate(d, d);
  }
  ZZ g = GCD(GCD(a, b), d);
  if (g != 1) {
    div(a, a, g);
    div(b, b, g);
    div(d, d, g);
  }

  rel_gen_a = to_long(a);
  rel_gen_b = to_long(b);
  rel_gen_d = to_long(d);
}

// Debug Tools
//  std::cout << "USING NUCOMP" << std::endl;
//  std::cout << "rgA is " << rgA << std::endl;
//  std::cout << "rgB is " << rgB << std::endl;
//  std::cout << "rgC is " << rgC << std::endl;
//  std::cout << "Ca is " << Ca << std::endl;
//  std::cout << "Cb is " << Cb << std::endl;
//  std::cout << "Cc is " << Cc << std::endl;
//  std::cout << "NUCOMP: RG1 is " << RelativeGenerator->conv_RR() << std::endl;
//  std::cout << "NUCOMP: RG2 is " <<
//  C.get_QO()->get_red_best()->get_RelativeGenerator()->conv_RR() << std::endl;
//  std::cout << "NUCOMP: RGf is " << RelativeGenerator->conv_RR() << std::endl;
//  std::cout << "NUCOMP: distance is " << log(RelativeGenerator->conv_RR()) <<
//  std::endl;
