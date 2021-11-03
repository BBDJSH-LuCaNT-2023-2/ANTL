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

  // normalize ideal
  if (!A.is_normal()) {
    normalize(A);
  }

  // reduce
  while (!A.is_reduced()) {
    a = A.get_a();
    b = A.get_b();
    c = A.get_c();

    A.assign(c, -1*b, a);
    normalize(A);
  }

  //account for special case
  if ((a == c) && (b < 0)) {
    NTL::negate(b,b);
    A.set_b(b);
  }
}

template <> void ReducePlainReal<ZZ>::normalize(QuadraticIdealBase<ZZ> & A) {
  static ZZ a, b, c, a2, delta, rootDelta, temp, s;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // delta = b^2 - 4ac
  mul(temp, a, c);
  mul(temp, temp, 4);
  sqr(delta, b);
  sub(delta, delta, temp);

  rootDelta = SqrRoot(delta);

  if(a <= rootDelta) {
    mul(a2, 2, abs(a));

    // Computing s, the normalizing integer,  per [BV07, pg. 108]
    sub(temp, rootDelta, b);
    div(s, temp, a2);
    mul(s, s, sign(a));

    //c = a*s^2 + b*s + c
    mul(temp, s, s);
    mul(temp, temp, a);
    add(c, c, temp);
    mul(temp, b, s);
    add(c, c, temp);

    //b = b + 2sa
    mul(temp, a, 2);
    mul(temp, temp, s);
    add(b, b, temp);

    A.assign(a, b, c);
  }

  else {
    mul(a2, 2, abs(a));

    // Computing s, the normalizing integer,  per [BV07, pg. 108]
    sub(temp, abs(a), b);
    div(s, temp, a2);
    mul(s, s, sign(a));

    //c = a*s^2 + b*s + c
    mul(temp, s, s);
    mul(temp, temp, a);
    add(c, c, temp);
    mul(temp, b, s);
    add(c, c, temp);

    //b = b + 2sa
    mul(temp, a, 2);
    mul(temp, temp, s);
    add(b, b, temp);

    A.assign(a, b, c);
  }
}
