#include <ANTL/Quadratic/Multiply/MultiplyNucomp_Opt.hpp>
#include <NTL/RR.h>

template <>
void MultiplyNucompOpt<ZZ>::construct_relative_generator(
    ZZ &rel_gen_a, ZZ &rel_gen_b, ZZ &rel_gen_d, QuadraticIdealBase<ZZ> &C,
    ZZ OB, ZZ BB, ZZ S);

template <>
void MultiplyNucompOpt<ZZ>::init(const ZZ &delta_in, const ZZ &h_in,
                                 long g_in) {
  init(delta_in, h_in, 0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}

template <>
void MultiplyNucompOpt<ZZ>::multiply(QuadraticIdealBase<ZZ> &C,
                                     const QuadraticIdealBase<ZZ> &A,
                                     const QuadraticIdealBase<ZZ> &B) {

  Delta = C.get_QO()->get_discriminant();
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));

  static ZZ a1, a2, b1, b2, c2, Ca, Cb, Cc, ss, m;
  static ZZ SP, S, v1, u2, v2, K, T, temp;
  static ZZ R1, R2, C1, C2, M1, M2;

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
  add(ss, b1, b2);
  RightShift(ss, ss, 1);

  sub(m, b1, b2);
  RightShift(m, m, 1);

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT(SP, v1, a2, a1);

  // K = v1 (b1 - b2) / 2 (mod L)
  mul(K, m, v1);
  rem(K, K, a1);

  S = 1;
  if (!IsOne(SP)) {
    XGCD(S, u2, v2, SP, ss);

    // K = u2 K - v2 c2 (mod L)
    mul(K, K, u2);
    mul(temp, v2, c2);
    sub(K, K, temp);

    if (!IsOne(S)) {
      div(a1, a1, S);
      div(a2, a2, S);
      mul(c2, c2, S);
    }

    rem(K, K, a1);
  }

  // N = a2;  L = a1;

  // check if NUCOMP steps are required
  if (a1 <= NC_BOUND && false) {
    // compute with regular multiplication formula (result will be reduced)

    // T = NK
    mul(T, a2, K);

    // C.a = A.a B.a / d^2 = NL
    mul(Ca, a2, a1);

    // C.b = b2 + 2 a2 K = b2 + 2 T
    LeftShift(Cb, T, 1);
    add(Cb, Cb, b2);

    // C.c = (S c2 + K (b2 + T)) / L;
    add(Cc, b2, T);
    mul(Cc, Cc, K);
    add(Cc, Cc, c2);
    div(Cc, Cc, a1);

    // Set a, b, c (DO NOT ASSIGN/NORMALIZE)
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
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // M1 = (N R1 + (b1 - b2) C1 / 2) / L  (T = N R1)
    mul(T, a2, R1);
    mul(M1, m, C1);
    add(M1, M1, T);
    div(M1, M1, a1);

    // M2 = (R1(b1 + b2)/2 - c2 S C1) / L
    mul(M2, ss, R1);
    mul(temp, c2, C1);
    sub(M2, M2, temp);
    div(M2, M2, a1);

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mul(Ca, R1, M1);
    mul(temp, C1, M2);
    if (C1 > 0)
      sub(Ca, Ca, temp);
    else
      sub(Ca, temp, Ca);

    // C.b = 2 (N R1 - C.a C2) / C1 - b1
    Cb = (a2 * R1 + C2 * Ca) << 1;
    Cb /= C1;
    Cb -= b2;

    // C.c = (C.b^2 - Delta) / 4 C.a
    sqr(Cc, Cb);
    sub(Cc, Cc, Delta);
    div(Cc, Cc, Ca);
    RightShift(Cc, Cc, 2);


    // Set a, b, c (DO NOT ASSIGN/NORMALIZE)
    C.set_a(Ca);
    C.set_b(Cb);
    C.set_c(Cc);

    // Partial update of the distance
    // arb_add(qie->distance, qie->distance, qie->distance, precision);
  }

  // Reduce and get the coefficients of the relative generator
  C.reduce();

  // rrelative_generator
  ZZ rel_gen_a, rel_gen_b, rel_gen_d;
  RR relative_generator;

  construct_relative_generator(rel_gen_a, rel_gen_b, rel_gen_d, C, abs(C2),
                               abs(C1), S);

  RelativeGenerator->set_abd(rel_gen_a, rel_gen_b, rel_gen_d);
  RelativeGenerator->invert();
  if (RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, ZZ(-1));
  }

}

//
// construct relative generator
// computes the coefficients a, b, d of the relative generator
// gamma = (a + b*sqrt(Delta))/d such that A^2 = (gamma) B
// Pass OB=1 and BB=0
//
template <>
void MultiplyNucompOpt<ZZ>::construct_relative_generator(
    ZZ &rel_gen_a, ZZ &rel_gen_b, ZZ &rel_gen_d, QuadraticIdealBase<ZZ> &C,
    ZZ OB, ZZ BB, ZZ S) {
  static ZZ NB;
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
