/**
 * @file qo_reduce_fast_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_reduce_fast class for GF2EX. 
 */

#include <ANTL/Quadratic/Reduce/ReduceFast.hpp>

template <> void ReduceFast<GF2EX>::reduce (QuadraticIdealBase<GF2EX> & A) {
  GF2EX a, b, c, q, r, temp;

  a = A.get_a();
  b = A.get_b();
  mul(c,A.get_c(),LeadCoeff (a));
  MakeMonic (a);

  if (deg(a) <= genus) {
    // already reduced - just normalize the ideal if necessary
    if (deg (b) >= deg (a))
      {
        // q = b/a
        DivRem (q, r, b, a);

        // c -= q(b + r + hx)
        add(temp,b,r);
        add(temp,temp,hx);
        mul(temp,temp,q);
        add(c,c,temp);

        // b = b % a
        b = r;
      }
  }
  else {
    // reduce
    GF2EX RR,R,BB,B,nB,oa;

    // run partial XGCD, using bound (deg(a) + g) / 2
    long N = (deg(a) + genus + 1) >> 1;
    bool even = !((deg(a) - genus) & 1);

    RR = b;
    R = a;
    oa = a;

    XGCD_PARTIAL_REDUCE(RR,R,BB,B,N,even);

    // a = (R^2 + R B h + Delta B^2) / oa
    mul(a,R,hx);
    mul(temp,Delta,B);
    add(a,a,temp);
    //    mul(a,a,B);
    //   sqr(temp,R);
    ::MulExact(a,a,B,deg(oa));
    ::SqrExact(temp,R,deg(oa));
    add(a,a,temp);
    div(a,a,oa);

    // b = h + (R + a BB) / B
    //    mul(temp,a,BB);
    ::MulExact(temp,a,BB,deg(B));
    add(b,R,temp);
    div(b,b,B);
    add(b,b,hx);

    // normalize and compute c
    MakeMonic(a);
    rem(b,b,a);

    //    sqr(c,b);
    //    mul(temp,b,hx);
    ::SqrExact(c,b,deg(a));
    ::MulExact(temp,b,hx,deg(a));
    add(c,c,temp);
    add(c,c,Delta);
    div(c,c,a);
  }

  A.assign(a,b,c);

#ifdef DEBUG_FASTRED
  if (deg(a) > genus) {
    cout << "ERROR (qo_reduce_fast):  not reduced" << endl;
    cout << A << endl;
    exit(1);
  }
#endif
}
