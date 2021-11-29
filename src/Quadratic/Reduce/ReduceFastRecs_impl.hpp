/**
 * @file qo_reduce_fast_recs_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of the qo_reduce_fast_recs class 
 * (for odd char base fields).
 */

template < class T > 
void 
qo_reduce_fast_recs<T>::reduce(QuadraticIdealBase<T> & A)
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
    T RR,R,AA,A,nA,BB,B,nB,E;
    bool flag;

    // run partial XGCD, using bound (deg(a) - g) / 2
    long N = (deg(a) - genus) >> 1;
    bool even = !((deg(a) - genus) & 1);

    clear(AA);
    set(A);

    set (BB);
    clear (B);

    RR = b;
    R = a;

    flag = true;

    // RR = q*R + r
    DivRem (q, r, RR, R);

    // nB = BB + q*B
    mul(nB,B,q);
    add(nB,nB,BB);

    while (!((deg(nB) == N) && even) && (deg(nB) <= N)) {
      BB = B;
      B = nB;

      mul(nA,A,q);
      add(nA,nA,AA);
      AA = A;
      A = nA;

      RR = R;
      R = r;

      flag = !flag;

      // RR = q*R + r
      DivRem(q,r,RR,R);

      // nB = BB + q*B
      mul(nB,B,q);
      add(nB,nB,BB);
    }


    // compute E_i = c_0 B_i - b_0 A_i
    mul(E,c,B);
    mul(temp,b,A);
    sub(E,E,temp);

    if (flag) {
      // a = A R + E B
      mul(a,A,R);
      mul(temp,E,B);
      add(a,a,temp);

      // b = -(AA R + E BB)
      mul(b,AA,R);
      mul(temp,E,BB);
      add(b,b,temp);
      NTL::negate(b,b);
    }
    else {
      // a = A R - E B
      mul(a,A,R);
      mul(temp,E,B);
      sub(a,a,temp);

      // b = E BB - AA R
      mul(b,E,BB);
      mul(temp,AA,R);
      sub(b,b,temp);
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
    cout << "ERROR (qo_reduce_fast_recs):  not reduced" << endl;
    cout << A << endl;
    exit(1);
  }
#endif
}
