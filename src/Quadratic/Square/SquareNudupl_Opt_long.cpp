#include <ANTL/Quadratic/Square/SquareNudupl_Opt.hpp>
#include <NTL/RR.h>

template <>
void SquareNuduplOpt<long>::construct_relative_generator(long &rel_gen_a, long &rel_gen_b, long &rel_gen_d,
                                  QuadraticIdealBase<long> &C, long OB, long BB,
                                  long S);

template <>
void SquareNuduplOpt<long>::init(const long &delta_in, const long &h_in, long g_in) {
  init(delta_in, h_in, 0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}

template <>
void SquareNuduplOpt<long>::square(QuadraticIdealBase<long> &C,
                                 const QuadraticIdealBase<long> &A) {
  // temporary computation of delta and nc_bound
  Delta = C.get_QO()->get_discriminant();
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));

  static long a1, b1, c1, Ca, Cb, Cc;
  static long S, v1, K, T;
  static long R1, R2, C1, C2, M2, temp;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();

  // solve S = v1 b1 + u1 a1 (only need v1)
  XGCD_LEFT(S, v1, b1, a1);

  // K = -v1 c1 (mod L)
  K = -(v1 * c1);

  if (S != 1) {
    a1 /= S;
    c1 *= S;
  }

  K %= a1;
  if (K < 0)
    K += a1;

  // N = L = a1
  // check if NUCOMP steps are required
  if (a1 <= NC_BOUND) {

    // compute with regular squaring formula (result will be reduced)

    // T = NK
    T = a1 * K;

    // C.a = a1^2 / S^2 = N^2
    // C.b = b1 + 2 a1 K = b1 + 2 T
    // C.c = (S c1 + K (b1 + T)) / L;
    Ca = a1 * a1;
    Cb = b1 + (T << 1);
    Cc = (c1 + (K * (b1 + T))) / a1;

    // Assign and normalize
    C.assign(Ca, Cb, Cc);

    // Partial update of the distance
    // arb_add(qie->distance, qie->distance, qie->distance, precision);

    C2 = 1;
    C1 = 0;
  } else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2 = a1;
    R1 = K;
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // M1 = R1

    // M2 = (R1 b1 - c1 S C1) / L
    M2 = (R1 * b1 - c1 * C1) / a1;

    // C.a = (-1)^(i-1) (R1^2 - C1 M2)
    Ca = R1 * R1;
    temp = C1 * M2;
    if (C1 > 0)
      Ca -= temp;
    else
      Ca = temp - Ca;

    // C.b = 2 (N R1 - C.a C2) / C1 - b1
    Cb = (a1 * R1 + C2 * Ca) << 1;
    Cb /= C1;
    Cb -= b1;

    // C.c = (C.b^2 - Delta) / 4 C.a
    Cc = (Cb * Cb - Delta) / Ca;
    Cc >>= 2;

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
  long rel_gen_a, rel_gen_b, rel_gen_d;
  RR relative_generator;

  construct_relative_generator(rel_gen_a, rel_gen_b, rel_gen_d, C, abs(C2),
                               abs(C1), S);

  RelativeGenerator->set_abd(rel_gen_a, rel_gen_b, rel_gen_d);
  RelativeGenerator->invert();
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
void SquareNuduplOpt<long>::construct_relative_generator(long &rel_gen_a, long &rel_gen_b, long &rel_gen_d,
                                  QuadraticIdealBase<long> &C, long OB, long BB,
                                  long S) {
  static long NB;
  bool con_rel_gen_dbg = false;

  if(con_rel_gen_dbg) {
    // printf("\n-->--> construct_relative_generator:\n");
    // printf("OB=%ld, BB=%ld\n",OB,BB);
    std::cout << "\n-->--> construct_relative_generator:" << std::endl;
    std::cout << "OB=" << OB << ", BB=" << BB << std::endl;
  }


  for (long i = 0; i < C.get_num_q(); ++i) {
    NB = C.get_qlist_i(i) * BB + OB;
    OB = BB;
    BB = NB;

    if(con_rel_gen_dbg) {
      // printf("&&& q=%ld, B=%ld\n",qie->qlist[i],BB);
      std::cout << "&&& q=" << C.get_qlist_i(i) << ", B=" << BB << std::endl;
    }
  }

  rel_gen_a = S * ((OB * C.get_a() << 1) + BB * C.get_b());
  rel_gen_b = -S * BB;
  rel_gen_d = (C.get_a() << 1);

  if(con_rel_gen_dbg) {
    // printf("OB=%ld, BB=%ld, a=%ld, b=%ld, S=%ld\n",OB,BB,qie->a,qie->b,S);
    // printf("relgen_a=%ld, relgen_b=%ld, relgen_d=%ld\n",*rel_gen_a,*rel_gen_b,*rel_gen_d);
    // printf("-->--> done construct_relative_generator\n\n");

    std::cout << "OB=" << OB << ", BB=" << BB << ", a=" << C.get_a() << ", b=" << C.get_b() << ", S=" << S<< std::endl;
    std::cout << "relgen_a=" << rel_gen_a << ", relgen_b=" << rel_gen_b << ", relgen_d=" << rel_gen_d << std::endl;
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
