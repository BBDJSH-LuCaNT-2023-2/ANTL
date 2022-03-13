/**
 * @file qo_nucube_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucube class for ZZ.
  */

#include <ANTL/Quadratic/Cube/CubeNucube.hpp>
#include <NTL/RR.h>


template <> void CubeNucube<ZZ>::init(const ZZ & delta_in, const ZZ & h_in, long g_in) {
  CubeStrategy<ZZ>::init(delta_in,h_in,0);
  sqrt_delta = FloorToZZ(sqrt(abs(to_RR(delta_in))));
}


template <> void CubeNucube<ZZ>::cube(QuadraticIdealBase<ZZ> &C, const QuadraticIdealBase<ZZ> &A) {
  Delta = C.get_QO()->get_discriminant();
  sqrt_delta = FloorToZZ(sqrt(abs(to_RR(Delta))));

  static ZZ a, b, c, Ca, Cb, Cc;
  static ZZ SP, S, v1, u2, v2, N, K, L, T, temp, temp2;
  static ZZ B, R1, R2, C1 = ZZ(-1), C2 = ZZ(0), BB, M1, M2, B1, B2, rgA, rgB, rgC;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // solve SP = v1 b + u1 a (only need v1)
  XGCD_LEFT (SP, v1, b, a);

  std::cout << "NUCUBE: v1 is " << v1 << std::endl;

  if (IsOne(SP)) {
    // N = a, S = 1, L = a^2
    N = a;
    S = 1;

    sqr(L,a);

    // K = c v1 (v1(b - a c v1) - 2) mod L
    rem(v1,v1,L);
    mul(temp,v1,c);
    rem(temp,temp,L);

    MulMod(K,temp,a,L);
    SubMod(K,b,K,L);
    MulMod(K,K,v1,L);
    SubMod(K,K,2,L);
    MulMod(K,K,temp,L);
  }
  else {
    // S = u2 (a SP) + v2 (b^2 - ac)
    mul(SP,SP,a);

    sqr(temp,b);
    mul(T,a,c);
    sub(temp,temp,T);

    XGCD(S,u2,v2,SP,temp);

    // N = a/S
    div(N,a,S);

    // L = N a
    mul(L,N,a);

    // K = -c(v1 u2 a + v2 b) mod L
    rem(v2,v2,L);
    mul(K,v1,u2);
    rem(K,K,L);
    MulMod(K,K,a,L);
    MulMod(temp,v2,b,L);
    AddMod(K,K,temp,L);
    MulMod(K,K,c,L);
    sub(K,L,K);

    // C = Sc
    mul(c,c,S);
  }

  // Compute NUCOMP termination bound
  //   mul(B, sqrt_delta, a);
  //   RightShift(B,B,1);
  //   SqrRoot(B, B);

    ZZ NC_BOUND = FloorToZZ(sqrt(to_RR(a)) * (sqrt(sqrt(abs(to_RR(Delta)/4)))));

  if (L < NC_BOUND) {
    std::cout << "NUCUBE: PLAIN" << std::endl;
    // compute with regular cubing formula (result will be reduced)

    // T = NK
    mul(T,N,K);

    // C.a = N L
    mul(Ca,N,L);

    // C.b = b + 2 T
    LeftShift(Cb,T,1);
    add(Cb,Cb,b);

    // C.c = (S c + K (T + b)) / L
    add(Cc,T,b);
    mul(Cc,Cc,K);
    add(Cc,Cc,c);
    div(Cc,Cc,L);
  }
  else {
    std::cout << "NUCUBE: NUCOMP" << std::endl;
    // use NUCOMP formulas

    // Execute partial reduction
    R2=L; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // T = N K
    MulMod(T,N,K,L);

    // M1 = (N R1 + T C1) / L  (temp = N R1)
    mul(temp,N,R1);
    mul(M1,T,C1);
    add(M1,M1,temp);
    div(M1,M1,L);

    // M2 = (R1(b + T) - c S C1) / L
    add(M2,b,T);
    mul(M2,M2,R1);
    mul(temp2,c,C1);
    sub(M2,M2,temp2);
    div(M2,M2,L);

    std::cout << "NUCUBE: R1 is " << R1 << std::endl;
    std::cout << "NUCUBE: M1 is " << M1 << std::endl;

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mul(Ca,R1,M1);
    mul(temp2,C1,M2);
    std::cout << "NUCUBE: R1M1 is " << Ca << std::endl;
    std::cout << "NUCUBE: C1M2 is " << temp2 << std::endl;
    if (sign(C1) < 0)
      sub(Ca,Ca,temp2);
    else
      sub(Ca,temp2,Ca);

    // C.b = 2 (N R1 - C.a C2) / C1 - b (mod 2a)
    mul(Cb,Ca,C2);
    sub(Cb,temp,Cb);
    LeftShift(Cb,Cb,1);
    div(Cb,Cb,C1);
    sub(Cb,Cb,b);
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

  std::cout << "Ca is " << Ca << std::endl;
  std::cout << "Cb is " << Cb << std::endl;
  std::cout << "Cc is " << Cc << std::endl;

  std::cout << "C1 is " << C1 << std::endl;
  std::cout << "C2 is " << C2 << std::endl;

  std::cout << "B1 is " << B1 << std::endl;
  std::cout << "B2 is " << B2 << std::endl;

  std::cout << "S is " << S << std::endl;

  rgA = S*(2*Ca*B1 + B2*Cb);
  rgB = -S*B2;
  rgC = 2*Ca;

  std::cout << "rgA is " << rgA << std::endl;
  std::cout << "rgB is " << rgB << std::endl;
  std::cout << "rgC is " << rgC << std::endl;

  RelativeGenerator->set_abd(rgA, rgB, rgC);
  RelativeGenerator->invert();
  if(RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, ZZ(-1));
  }

  std::cout << "Ca is " << Ca << std::endl;
  std::cout << "Cb is " << Cb << std::endl;
  std::cout << "Cc is " << Cc << std::endl;

  C.assign(Ca,Cb,Cc);
  C.reduce();

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
