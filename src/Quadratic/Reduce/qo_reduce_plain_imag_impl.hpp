/**
 * @file qo_reduce_plain_impl.hpp
 * @author Michael Jacobson
 * @remark Generic implementation of basic ideal reduction (odd char base types).
 */

//
// reduce
//
// Task:
//      reduces the ideal
//

template < class T > 
void 
qo_reduce_plain_imag<T>::reduce(QuadraticIdealBase<T> & A)
{
  T a, b, c, na, nb, q, r, a2, temp;

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
    // reduce
    while (deg (a) > genus)
      {
	na = c;

	// -b = q*na + nb
	DivRem (q, nb, -b, c);

	// c = a - q*(nb - b) = a + q*(b - nb);
	sub(c,b,nb);
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
