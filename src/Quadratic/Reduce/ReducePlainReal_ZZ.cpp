/**
 * @file qo_reduce_plain_ZZ.cpp
 * @author Michael Jacobson
 * @remark Basic ideal reduction specialization (ZZ base type).
 */

#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

// reduce
//
// Task: reduces the ideal

template <> void ReducePlainReal<ZZ>::reduce(QuadraticIdealBase<ZZ> & A) {
  static ZZ a, b, c, na, nb, q, r, a2, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // normalize ideal
  if (b <= -a || b > a) {
    LeftShift(a2,a,1);
  
    // q = b/2a
    DivRem(q, r, b, a2);

    if (r > a) {
      sub(r,r,a2);
      ++q;
    }

    // c -= q (b + r) / 2
    add(temp,b,r);
    RightShift(temp,temp,1);
    mul(temp,temp,q);
    sub(c,c,temp);

    // b = r
    b = r;
  }

  // reduce
  while (a > c) {
    na = c;

    LeftShift(a2,na,1);

    // -b = 2q * na + nb
    NTL::negate(temp,b);
    DivRem (q, nb, temp, a2);

    if (nb > na)
      {
    sub(nb,nb,a2);
    ++q;
      }

    // c = a - q * (nb - b)/2
    sub(temp,nb,b);
    RightShift(temp,temp,1);
    mul(temp,temp,q);
    sub(c,a,temp);

    b = nb;
    a = na;
  }

  // account for special case
  if ((a == c) && (b < 0))
    NTL::negate(b,b);

  A.assign(a,b,c);
}
