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
  static long a, b, c;

  // normalize ideal
  if (!A.is_normal()) {
    A.normalize();
  }

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  // reduce
  while (!A.is_reduced()) {
    A.assign(c, -b, a);
    A.normalize();

    a = A.get_a();
    b = A.get_b();
    c = A.get_c();
  }

  if (a < 0) {
    A.set_a(-a);
    A.set_c(-c);
  }

  //account for special case
  if ((a == c) && (b < 0)) {
    b = -b;
    A.set_b(b);
  }
}
