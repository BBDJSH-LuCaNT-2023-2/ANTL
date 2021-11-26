/**
 * @file qo_reduce_plain_ZZ.cpp
 * @author Michael Jacobson
 * @remark Basic ideal reduction specialization (ZZ base type).
 */

#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

// reduce
//
// Task:
//      reduces the ideal

template <>
void ReducePlainReal<ZZ>::reduce(QuadraticIdealBase<ZZ> & A) {
  static ZZ a, b, c;

  // normalize ideal
  if (!A.is_normal()) {
    A.normalize();
  }

  bool debug = false;
  // A variable Xn refers to X_{i+n} E.g. R2 refers to R_{i+2}
  ZZ R, Q, P, B;
  ZZ q, r;

  RR RelativeGenerator;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  while() {
    DivRem(q, r, b + FloorRootDelta, 2*a);

    if(debug) {
      std::cout << "a is " << a << std::endl;
      std::cout << "b is " << b << std::endl;
      std::cout << "c is " << c << std::endl;
      std::cout << "q is " << q << std::endl;
      std::cout << "r is " << r << std::endl;
    }

    R = -a;
    P = FloorRootDelta - r;
    Q = q*((b - P)/2) - c;

    if(debug) {
    std::cout << "Q is " << Q << std::endl;
    std::cout << "P is " << P << std::endl;
    std::cout << "R is " << R << std::endl;
  }

  }
  RelativeGenerator = inv(abs((to_RR(P) - sqrt(to_RR(Delta))) / to_RR(2*Q)));

  if(debug) {
    std::cout << "RG is " << RelativeGenerator << std::endl;
    std::cout << "distance is " << log(RelativeGenerator) << std::endl;
  }


  A.assign(Q, P, R);
  distance = log(RelativeGenerator);
}
