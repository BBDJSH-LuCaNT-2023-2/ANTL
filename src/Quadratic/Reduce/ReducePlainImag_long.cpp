/**
 * @file qo_reduce_plain_long.cpp
 * @author Michael Jacobson
 * @remark Basic ideal reduction specialization (long base type).
 */

#include <ANTL/Quadratic/Reduce/ReducePlainImag.hpp>

// reduce
//
// Task: reduces the ideal

template <> void ReducePlainImag<long>::reduce(QuadraticIdealBase<long> & A) {
  static long a, b, c, na, nb, q, r, a2, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // normalize ideal
  if (b <= -a || b > a) {
    a2 = a << 1;  
  
    // q = b/2a
    DivRem(q, r, b, a2);

    if (r > a) {
      r -= a2;
      ++q;
    }

    // c -= q (b + r) / 2
    c -= q*(b+2) >> 1;

    // b = r
    b = r;
  }

  // reduce
  while (a > c) {
    na = c;

    a2 = na << 1;

    // -b = 2q * na + nb
    DivRem (q, nb, -b, a2);

    if (nb > na)
      {
	nb -= a2;
	++q;
      }

    // c = a - q * (nb - b)/2
    c = a - q*(nb - b) >> 1;

    b = nb;
    a = na;
  }

  // account for special case
  if ((a == c) && (b < 0))
    b = -b;

  A.assign(a,b,c);
}
