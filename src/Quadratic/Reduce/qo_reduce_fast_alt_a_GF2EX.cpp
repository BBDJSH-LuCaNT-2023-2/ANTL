/**
 * @file qo_reduce_fast_alt_a_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_reduce_fast_alt_a class for GF2EX. 
 */

#include <Quadratic/Reduce/qo_reduce_fast_alt_a.hpp>

template <> 
void 
qo_reduce_fast_alt_a<GF2EX>::reduce (QuadraticIdealBase<GF2EX> & A)
{
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
    GF2EX RR,R,BB,B,nB,A,E;

    // run partial XGCD, using bound (deg(a) - g) / 2
    long N = (deg(a) + genus + 1) >> 1;
    bool even = !((deg(a) - genus) & 1);

    RR = b;
    R = a;

    XGCD_PARTIAL_REDUCE(RR,R,BB,B,N,even);

    // A = (B b0 + R) / a0
    //    mul(A,B,b);
    ::MulExact(A,B,b,deg(a));
    add(A,A,R);
    div(A,A,a);

    // E = c0 B + (h + b0) A
    mul(E,c,B);
    add(temp,hx,b);
    mul(temp,temp,A);
    add(E,E,temp);

    // a = A R + E B
    mul(a,A,R);
    mul(temp,E,B);
    add(a,a,temp);

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
    cout << "ERROR (qo_reduce_fast_alt_a):  not reduced" << endl;
    cout << A << endl;
    exit(1);
  }
#endif
}
