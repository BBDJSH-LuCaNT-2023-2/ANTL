/**
 * @file qo_nucube_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_nucube class (for odd char base fields).
 */

template <class T> void CubeNucube<T>::cube(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A) {
  T a, b, c, b2, Ca, Cb, Cc;
  T SP, S, v1, u2, v2, N, K, L, TT;
  T B, R1, R2, C1, C2, BB, M1, M2, temp, temp2;
  long BOUND;
  bool flag;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();
  add(b2,b,b);

  // solve SP = v1 (2 b) + u1 a (only need v1)
  XGCD_LEFT (SP, v1, b2, a);

  if (IsOne(SP)) {
    // N = a
    N = a;

    // L = a^2
    sqr(L,a);

    // K = c v1 (v1 (2 b - a c v1) - 2)  (mod L)
    mul(temp,v1,c);
    rem(temp,temp,L);
    mul(K,temp,a);
    rem(K,K,L);
    sub(K,b2,K);
    mul(K,K,v1);
    sub(K,K,2);
    rem(K,K,L);
    mul(K,K,temp);
    rem(K,K,L);
  }
  else {
    // S = u2 (a SP) + v2 (3 b^2 + D)
    mul(SP,SP,a);

    sqr(temp,b);
    mul(temp,temp,3);
    add(temp,temp,Delta);

    XGCD(S,u2,v2,SP,temp);

    // N = a/S
    div(N,a,S);

    // L = N a
    mul(L,N,a);

    // K = -c(v1 u2 a + 2 v2 b)  (mod L)
    mul(K,v1,u2);
    rem(K,K,L);
    mul(K,K,a);
    rem(K,K,L);
    mul(temp,v2,b2);
    rem(temp,temp,L);
    add(K,K,temp);
    mul(K,K,c);
    rem(K,K,L);

    // C = Sc
    mul(c,c,S);
  }

  // Compute NUCOMP termination bound
  BOUND = (deg(a) + genus) >> 1;

  // check if NUCOMP steps are required
  if (deg(L) <= BOUND) {
    // compute with regular multiplication formula (result will be reduced)

    // T = N K
    mul(TT,N,K);

    // C.a = N L
    mul(Ca,N,L);

    // C.b = b + T
    add(Cb,TT,b);

    // C.c = (S c + K (2 b + T)) / L
    add(Cc,TT,b2);
    //    mul(Cc,Cc,K);
    ::MulExact(Cc,Cc,K,deg(L));
    add(Cc,Cc,c);
    div(Cc,Cc,L);
  }
  else {
    // use NUCOMP formulas

    // Execute partial reduction
    R2=L; R1=K;
    XGCD_PARTIAL(R2, R1, C2, C1, BOUND, flag);

    // T = N K
    MulMod(TT,N,K,L);

    // M1 = (N R1 + TT C1) / L  (T = N R1)
    mul(temp,N,R1);
    //  mul(M1,TT,C1);
    ::MulExact(M1,TT,C1,deg(L));
    add(M1,M1,temp);
    div(M1,M1,L);

    // M2 = (R1(2 b + TT) - c2 S C1) / L
    add(M2,b2,TT);
    //   mul(M2,M2,R1);
    //    mul(temp2,c,C1);
    ::MulExact(M2,M2,R1,deg(L));
    ::MulExact(temp2,c,C1,deg(L));
    sub(M2,M2,temp2);
    div(M2,M2,L);

    // C.a = (-1)^(i-1) (R1 M1 - C1 M2)
    mul(Ca,R1,M1);
    mul(temp2,C1,M2);
    if (flag)
      sub(Ca,Ca,temp2);
    else
      sub(Ca,temp2,Ca);

    // C.b = (N R1 + C.a C2) / C1 - b (mod a)
    //    mul(C.b,C.a,C2);
    ::MulExact(Cb,Ca,C2,deg(C1));
    add(Cb,Cb,temp);
    div(Cb,Cb,C1);
    sub(Cb,Cb,b);
    rem(Cb,Cb,Ca);

    // C.c = (C.b^2 - Delta) / C.a
    //    sqr(Cc,Cb);
    ::SqrExact(Cc,Cb,deg(Ca));
    sub(Cc,Cc,Delta);
    div(Cc,Cc,Ca);
  }

  // normalize and reduce
  C.assign(a,b,c);
  C.reduce();
}
