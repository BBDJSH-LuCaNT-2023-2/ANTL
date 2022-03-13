/**
 * @file qo_nudupl_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nudupl class for ZZ.
 */

#include <ANTL/Quadratic/Square/SquareNudupl.hpp>
#include <NTL/RR.h>


template <> void SquareNudupl<ZZ>::init(const ZZ & delta_in, const ZZ & h_in, long g_in) {
  SquareStrategy<ZZ>::init(delta_in,h_in,0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}


template <> void SquareNudupl<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A) {
  Delta = C.get_QO()->get_discriminant();
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));

  static ZZ a1, b1, c1, Ca, Cb, Cc;
  static ZZ S, v1, K, T, temp;
  static ZZ R1, R2, C1 = ZZ(-1), C2 = ZZ(0), M2, B1, B2, rgA, rgB, rgC;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();

  // solve S = v1 b1 + u1 a1 (only need v1)
  XGCD_LEFT (S, v1, b1, a1);

  // K = -v1 c1 (mod L)
  mul(K,v1,c1);
  NTL::negate(K,K);

  if (!IsOne (S))
    {
      div(a1,a1,S);
      mul(c1,c1,S);
    }

  rem(K,K,a1);

  // N = L = a1;


  // check if NUCOMP steps are required
  if (a1 < NC_BOUND) {
    // compute with regular squaring formula (result will be reduced)

    // T = NK
    mul(T,a1,K);

    // C.a = A.a^2 / S^2 = N^2
    sqr(Ca,a1);

    // C.b = b1 + 2 a1 K = b1 + 2 T
    LeftShift(Cb,T,1);
    add(Cb,Cb,b1);

    // C.c = (S c1 + K (b1 + T)) / L;
    add(Cc,b1,T);
    mul(Cc,Cc,K);
    add(Cc,Cc,c1);
    div(Cc,Cc,a1);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // M1 = R1

    // M2 = (R1 b1 - c1 S C1) / L
    mul(M2,R1,b1);
    mul(temp,c1,C1);
    sub(M2,M2,temp);
    div(M2,M2,a1);

    // C.a = (-1)^(i-1) (R1^2 - C1 M2)
    sqr(Ca,R1);
    mul(temp,C1,M2);
    if (sign(C1) < 0)
      sub(Ca,Ca,temp);
    else
      sub(Ca,temp,Ca);

    // C.b = 2 (N R1 - C.a C2) / C1 - b1 (mod 2a)
    mul(Cb,Ca,C2);
    mul(temp,a1,R1);
    sub(Cb,temp,Cb);
    LeftShift(Cb,Cb,1);
    div(Cb,Cb,C1);
    sub(Cb,Cb,b1);
    rem(Cb,Cb,Ca << 1);

    // C.c = (C.b^2 - Delta) / 4 C.a
    sqr(Cc,Cb);
    sub(Cc,Cc,Delta);
    div(Cc,Cc,Ca);
    RightShift(Cc,Cc,2);

    if (Ca < 0) {
      NTL::negate(Ca,Ca);
      NTL::negate(Cc,Cc);
    }
  }

  B1 = ZZ(sign(Ca))*abs(C1);
  B2 = abs(C2);

  rgA = S*(2*Ca*B1 + B2*Cb);
  rgB = -S*B2;
  rgC = 2*Ca;

  RelativeGenerator->set_abd(rgA, rgB, rgC);
  RelativeGenerator->invert();
  if(RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, ZZ(-1));
  }

  // normalize and reduce
  C.assign(Ca,Cb,Cc);
  C.reduce();

  std::cout << "NUCOMP: RG1 is " << RelativeGenerator->conv_RR() << std::endl;
  std::cout << "NUCOMP: RG2 is " << C.get_QO()->get_red_best()->get_RelativeGenerator()->conv_RR() << std::endl;

  ANTL::mul(*RelativeGenerator, *RelativeGenerator, *C.get_QO()->get_red_best()->get_RelativeGenerator());
}
//Debug Tools
// std::cout << "USING NUCOMP" << std::endl;
// std::cout << "rgA is " << rgA << std::endl;
// std::cout << "rgB is " << rgB << std::endl;
// std::cout << "rgC is " << rgC << std::endl;
// std::cout << "Ca is " << Ca << std::endl;
// std::cout << "Cb is " << Cb << std::endl;
// std::cout << "Cc is " << Cc << std::endl;
// std::cout << "NUCOMP: RG1 is " << RelativeGenerator->conv_RR() << std::endl;
// std::cout << "NUCOMP: RG2 is " << C.get_QO()->get_red_best()->get_RelativeGenerator()->conv_RR() << std::endl;
// std::cout << "NUCOMP: RGf is " << RelativeGenerator->conv_RR() << std::endl;
// std::cout << "NUCOMP: distance is " << log(RelativeGenerator->conv_RR()) << std::endl;
