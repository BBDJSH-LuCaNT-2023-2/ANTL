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
  static ZZ a, b, c;
  std::cout << "begin reduce!" << std::endl;
  std::cout << "(a, b, c) = (" << A.get_a() << ", " << A.get_b() << ", " << A.get_c() << ")" << std::endl;
  // normalize ideal
  if (!A.is_normal()) {
    A.normalize();
    std::cout << "(a, b, c) = (" << A.get_a() << ", " << A.get_b() << ", " << A.get_c() << ")" << std::endl;
  // normalize ideal
  }

  // reduce
  while (!A.is_reduced()) {
    a = A.get_a();
    b = A.get_b();
    c = A.get_c();

    A.assign(c, -1*b, a);
    A.normalize();
    std::cout << "(a, b, c) = (" << A.get_a() << ", " << A.get_b() << ", " << A.get_c() << ")" << std::endl;
  // normalize ideal
  }

  //account for special case
  if ((a == c) && (b < 0)) {
    b = -b;
    A.set_b(b);
  }
}
