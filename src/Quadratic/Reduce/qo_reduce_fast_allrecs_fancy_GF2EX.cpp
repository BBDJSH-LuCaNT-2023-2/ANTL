/**
 * @file qo_reduce_fast_allrecs_fancyGF2EX.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_reduce_fast_allrecs_fancy class for GF2EX. 
 */

#include <Quadratic/Reduce/qo_reduce_fast_allrecs_fancy.hpp>

template <> 
void 
qo_reduce_fast_allrecs_fancy<GF2EX>::reduce (QuadraticIdealBase<GF2EX> & A)
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
    GF2EX RR,R,AA,A,nA,BB,B,nB,EE,E,nE;

    // run partial XGCD, using bound (deg(a) - g) / 2
    long N = (deg(a) - genus) >> 1;
    bool even = !((deg(a) - genus) & 1);

    clear(AA);
    set(A);

    set (BB);
    clear (B);

    EE = c;
    add(E,hx,b);

    RR = b;
    R = a;

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

      mul(nE,E,q);
      add(nE,nE,EE);
      EE = E;
      E = nE;

      RR = R;
      R = r;

      // RR = q*R + r
      DivRem(q,r,RR,R);

      // nB = BB + q*B
      mul(nB,B,q);
      add(nB,nB,BB);
    }

    // a = A R + E B
    mul(a,A,R);
    mul(temp,E,B);
    add(a,a,temp);

    // b = h + AA R + E BB
    mul(b,AA,R);
    mul(temp,E,BB);
    add(b,b,temp);
    add(b,b,hx);

    // c = AA RR + EE BB
    mul(c,AA,RR);
    mul(temp,EE,BB);
    add(c,c,temp);

    // normalize
    c *= LeadCoeff(a);
    MakeMonic(a);

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

  A.assign(a,b,c);

#ifdef DEBUG_FASTRED
  if (deg(a) > genus) {
    cout << "ERROR (qo_reduce_allrecs_fancy):  not reduced" << endl;
    cout << A << endl;
    exit(1);
  }
#endif
}
