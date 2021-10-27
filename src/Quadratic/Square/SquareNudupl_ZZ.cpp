/**
 * @file qo_nudupl_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nudupl class for ZZ.
 */

#include <ANTL/Quadratic/Square/SquareNudupl.hpp>
#include <NTL/RR.h>


template <> void SquareNudupl<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin) {
  SquareStrategy<ZZ>::init(Din,hin,0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}


template <> void SquareNudupl<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A) {
  static ZZ a1, b1, c1, Ca, Cb, Cc;
  static ZZ S, v1, K, T, temp;
  static ZZ R1, R2, C1, C2, M2;

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

  // normalize and reduce
  C.assign(Ca,Cb,Cc);
  //C.reduce();
}
