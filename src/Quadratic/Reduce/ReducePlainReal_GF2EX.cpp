/**
 * @file qo_reduce_plain_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Basic ideal reduction specialization (GF2EX base type).
 */

#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

// reduce_plain
//
// Task: reduces the ideal

template <> void ReducePlainReal<GF2EX>::reduce (QuadraticIdealBase<GF2EX> & A) {
  GF2EX a, b, c, na, nb, q, r, a2, temp;

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
    while (deg (a) > genus)
      {
	na = c;

	// b + hx = q*na + nb
	add(temp,b,hx);
	DivRem (q, nb, temp, c);

	// c = oa + q*(nb + b);
	add(c,nb,b);
	mul(c,c,q);
	add(c,c,a);

	b = nb;
	a = na;
      }

    // make sure a is monic
    mul(c,c,LeadCoeff(a));
    MakeMonic (a);
  }

  A.assign(a,b,c);
}
