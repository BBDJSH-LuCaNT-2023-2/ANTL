/**
 * @file qo_nudupl_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_nudupl class (for odd char base fields).
 */

template < class T > 
void
qo_nudupl<T>::square(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A)
{
  T a1, b1, c1, b2, Ca, Cb, Cc;
  T S, v1, K, TT;
  T R1, R2, C1, C2, M2, temp;
  long BOUND;
  bool flag;

  a1 = A.get_a();
  b1 = A.get_b();
  c1 = A.get_c();
  add(b2,b1,b1);

  // solve S = v1 (2 b1) + u1 a1 (only need v1)
  XGCD_LEFT (S, v1, b2, a1);

  // K = -v1 c1 (mod L)
  mul(K,v1,c1);
  NTL::negate(K,K);

  if (!IsOne (S))
    {
      div(a1,a1,S);
      mul(c1,c1,S);
    }

  rem(K,K,a1);

  // N = L = a1

  // check if NUCOMP steps are required
  if ((deg(a1) << 1) <= genus) {
    // compute with regular multiplication formula (result will be reduced)

    // T = NK
    mul(TT,a1,K);

    // C.a = A.a^2 / S^2 = N^2
    sqr(Ca,a1);

    // C.b = b1 + N K = b1 + T
    add(Cb,TT,b1);

    // C.c = (S c1 + K (2 b1 + T)) / L;
    add(Cc,b2,TT);
    //    mul(Cc,Cc,K);
    ::MulExact(Cc,Cc,K,deg(a1));
    add(Cc,Cc,c1);
    div(Cc,Cc,a1);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    if (deg(b1) == genus + 1)
      BOUND = (genus + 1) >> 1;
    else
      BOUND = (genus) >> 1;

    R2=a1; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, BOUND, flag);

    // M1 = R1

    // M2 = (2 R1 b1 - c1 S C1) / L
    //   mul(M2,b2,R1);
    //    mul(temp,c1,C1);
    ::MulExact(M2,b2,R1,deg(a1));
    ::MulExact(temp,c1,C1,deg(a1));
    sub(M2,M2,temp);
    div(M2,M2,a1);

    // C.a = (-1)^(i-1) (R1^2 - C1 M2)
    sqr(Ca,R1);
    mul(temp,C1,M2);
    if (flag)
      sub(Ca,Ca,temp);
    else
      sub(Ca,temp,Ca);

    // C.b = (N R1 + C.a C2) / C1 - b1 (mod a)
    //    mul(C.b,C.a,C2);
    //    mul(temp,a1,R1);
    ::MulExact(Cb,Ca,C2,deg(C1));
    ::MulExact(temp,a1,R1,deg(C1));
    add(Cb,temp,Cb);
    div(Cb,Cb,C1);
    sub(Cb,Cb,b1);
    rem(Cb,Cb,Ca);

    // C.c = (C.b^2 - Delta) / C.a
    //    sqr(Cc,Cb);
    ::SqrExact(Cc,Cb,deg(Ca));
    sub(Cc,Cc,Delta);
    div(Cc,Cc,Ca);
  }

  // normalize and reduce
  C.assign(Ca,Cb,Cc);
  C.reduce();
}
