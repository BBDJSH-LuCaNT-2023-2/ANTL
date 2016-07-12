/**
 * @file qo_reduce_fast_alt_a_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_reduce_fast_alt_a class 
 * (for odd char base fields).
 */

template < class T > 
void 
qo_reduce_fast_alt_a<T>::reduce(QuadraticIdealBase<T> & A)
{
  T a,b,c,q,r,temp;

  a = A.get_a();
  b = A.get_b();
  mul(c,A.get_c(),LeadCoeff(a));
  MakeMonic (a);

  if (deg(a) <= genus) {
    // already reduced - just normalize the ideal if necessary
    if (deg (b) >= deg (a))
      {
        // q = b/a
        DivRem (q, r, b, a);

        // c -= q(b + r)
        add(temp,b,r);
        mul(temp,temp,q);
        sub(c,c,temp);

        // b = b % a
        b = r;
      }
  }
  else {
    T RR,R,A,BB,B,nB,E;
    bool flag;

    // run partial XGCD, using bound (deg(a) - g) / 2
    long N = (deg(a) + genus + 1) >> 1;
    bool even = !((deg(a) - genus) & 1);

    RR = b;
    R = a;

    XGCD_PARTIAL_REDUCE(RR,R,BB,B,N,flag,even);

    if (flag) {
      // A = (B b0 + R) / a0
      //      mul(A,B,b);
      ::MulExact(A,B,b,deg(a));
      add(A,A,R);
      div(A,A,a);

      // E = c0 B - b0 A
      mul(E,c,B);
      mul(temp,b,A);
      sub(E,E,temp);

      // a = A R + E B
      mul(a,A,R);
      mul(temp,E,B);
      add(a,a,temp);

      // b = (R - a BB) / B
      //      mul(temp,a,BB);
      ::MulExact(temp,a,BB,deg(B));
      sub(b,R,temp);
      div(b,b,B);
    }
    else {
      // A = (B b0 - R) / a0
      //      mul(A,B,b);
      ::MulExact(A,B,b,deg(a));
      sub(A,A,R);
      div(A,A,a);

      // E = c0 B - b0 A
      mul(E,c,B);
      mul(temp,b,A);
      sub(E,E,temp);

      // a = A R - E B
      mul(a,A,R);
      mul(temp,E,B);
      sub(a,a,temp);
 
      // b = -(R + a*BB) / B
      //      mul(b,a,BB);
      ::MulExact(b,a,BB,deg(B));
      add(b,b,R);
      div(b,b,B);
      NTL::negate(b,b);
    }

    // normalize and compute c
    MakeMonic(a);
    rem(b,b,a);

    //    sqr(c,b);
    ::SqrExact(c,b,deg(a));
    sub(c,c,Delta);
    div(c,c,a);
  }

  A.assign(a,b,c);

#ifdef DEBUG_FASTRED
  if (deg(a) > genus) {
    cout << "ERROR (qo_reduce_fast_alt_a):  not reduced" << endl;
    cout << A << endl;
    exit(1);
  }
#endif
}
