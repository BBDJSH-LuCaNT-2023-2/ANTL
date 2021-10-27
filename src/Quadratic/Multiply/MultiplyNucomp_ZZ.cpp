/**
 * @file qo_nucomp_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_nucomp class for ZZ.
 */

#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>
#include <NTL/RR.h>

template <> void MultiplyNucomp<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin) {
  MultiplyStrategy<ZZ>::init(Din,hin,0);
  NC_BOUND = FloorToZZ(sqrt(sqrt(abs(to_RR(Delta)))));
}

template <> void MultiplyNucomp<ZZ>::multiply(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A, const QuadraticIdealBase<ZZ> & B) {
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
  }
  else {
    a1 = B.get_a();
    a2 = A.get_a();
    b1 = B.get_b();
    b2 = A.get_b();
    c2 = A.get_c();
  }

  // s = (b1 + b2)/2, m = (b1 - b2)/2
  add(ss,b1,b2);
  RightShift(ss,ss,1);

  sub(m,b1,b2);
  RightShift(m,m,1);

  // solve SP = v1 a2 + u1 a1 (only need v1)
  XGCD_LEFT (SP, v1, a2, a1);

  // K = v1 (b1 - b2) / 2 (mod L)
  mul(K,m,v1);
  rem(K,K,a1);

  if (!IsOne (SP))
    {
      XGCD (S, u2, v2, SP, ss);

      // K = u2 K - v2 c2 (mod L)
      mul(K,K,u2);
      mul(temp,v2,c2);
      sub(K,K,temp);

      if (!IsOne (S))
	{
	  div(a1,a1,S);
	  div(a2,a2,S);
	  mul(c2,c2,S);
	}

      rem(K,K,a1);
    }

  // N = a2;  L = a1;

  // check if NUCOMP steps are required
  if (a1 < NC_BOUND) {
    // compute with regular multiplication formula (result will be reduced)

    // T = NK
    mul(T,a2,K);

    // C.a = A.a B.a / d^2 = NL
    mul(Ca,a2,a1);

    // C.b = b2 + 2 a2 K = b2 + 2 T
    LeftShift(Cb,T,1);
    add(Cb,Cb,b2);

    // C.c = (S c2 + K (b2 + T)) / L;
    add(Cc,b2,T);
    mul(Cc,Cc,K);
    add(Cc,Cc,c2);
    div(Cc,Cc,a1);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, NC_BOUND);

    // M1 = (N R1 + (b1 - b2) C1 / 2) / L  (T = N R1)
    mul(T,a2,R1);
    mul(M1,m,C1);
    add(M1,M1,T);
    div(M1,M1,a1);

    // M2 = (R1(b1 + b2)/2 - c2 S C1) / L
    mul(M2,ss,R1);
    mul(temp,c2,C1);
    sub(M2,M2,temp);
    div(M2,M2,a1);

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mul(Ca,R1,M1);
    mul(temp,C1,M2);
    if (sign(C1) < 0)
      sub(Ca,Ca,temp);
    else
      sub(Ca,temp,Ca);

    // C.b = 2 (N R1 - C.a C2) / C1 - b2 (mod 2a)
    mul(Cb,Ca,C2);
    sub(Cb,T,Cb);
    LeftShift(Cb,Cb,1);
    div(Cb,Cb,C1);
    sub(Cb,Cb,b2);
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
  C.reduce();
}
