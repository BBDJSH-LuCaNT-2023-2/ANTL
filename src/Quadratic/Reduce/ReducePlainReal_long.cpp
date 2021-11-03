/**
 * @file qo_reduce_plain_long.cpp
 * @author Michael Jacobson
 * @remark Basic ideal reduction specialization (long base type).
 */

#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

// reduce
//
// Task: reduces the ideal
template <> void ReducePlainReal<long>::reduce(QuadraticIdealBase<long> & A) {
  static long a, b, c, na, nb, q, r, a2, temp;

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
    b = -1*b;
    A.set_b(b);
  }
}

template <> void ReducePlainReal<long>::normalize(QuadraticIdealBase<long> & A) {
  static long a, b, c, a2, delta, rootDelta, s;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  delta = b*b - 4*a*c;

  rootDelta = SqrRoot(delta);

  if(a <= rootDelta) {
    a2 = 2*abs(a);

    // Computing s, the normalizing integer,  per [BV07, pg. 108]
    s = sign(a)*((rootDelta - b)/a2);

    //c = a*s^2 + b*s + c
    c = a*s*s + b*s + c;

    //b = b + 2sa
    b = b + 2*s*a;

    A.assign(a, b, c);
  }

  else {
    a2 = 2*abs(a);

    // Computing s, the normalizing integer,  per [BV07, pg. 108]
    s = sign(a)*((abs(a) - b)/a2);

    //c = a*s^2 + b*s + c
    c = a*s*s + b*s + c;

    //b = b + 2sa
    b = b + 2*s*a;

    A.assign(a, b, c);
  }
}
