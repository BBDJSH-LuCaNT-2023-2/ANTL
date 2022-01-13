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
  ZZ a0, b0, c0, a1, b1, c1, q, r;
  ZZ B0 = ZZ(1), B1 = ZZ(0), BTemp;

  //mul(a0, a0, b1);
  bool debug = false;
  set(*RelativeGenerator);

  if(A.is_reduced()) {
    return;
  }

  a0 = A.get_a();
  b0 = A.get_b();
  c0 = A.get_c();

  if(debug) {
    std::cout << "a is " << a0 << std::endl;
    std::cout << "b is " << b0 << std::endl;
    std::cout << "c is " << c0 << std::endl;
    }

  // normalize ideal
  if (!A.is_normal()) {
    if(debug) {
      std::cout << "not normal! " << std::endl;
    }
    A.normalize();
  }

  a0 = A.get_a();
  b0 = A.get_b();
  c0 = A.get_c();

  if(debug) {
  std::cout << "a is " << a0 << std::endl;
  std::cout << "b is " << b0 << std::endl;
  std::cout << "c is " << c0 << std::endl;
  }

  if(A.is_reduced()) {
    return;
  }

  DivRem(q, r, b0 + FloorRootDelta, 2*a0);

    if(debug) {
      std::cout << "FloorRootDelta is " << FloorRootDelta << std::endl;
      std::cout << "a is " << a0 << std::endl;
      std::cout << "b is " << b0 << std::endl;
      std::cout << "c is " << c0 << std::endl;
      std::cout << "q is " << q << std::endl;
      std::cout << "r is " << r << std::endl;
    }

    c1 = -a0;
    b1 = FloorRootDelta - r;
    a1 = q*((b0 - b1)/2) - c0;

    BTemp = B1;
    B1 = q*B1 + B0;
    B0 = B1;

    a0 = a1;
    b0 = b1;
    c0 = c1;

    if(debug) {
      std::cout << "a1 is " << a1 << std::endl;
      std::cout << "b1 is " << b1 << std::endl;
      std::cout << "c1 is " << c1 << std::endl;

      std::cout << "B0 is " << B0 << std::endl;
      std::cout << "B1 is " << B1 << std::endl;
    }

  while(!(abs(FloorRootDelta - 2*a0) <= b0 - 1 && b0 <= FloorRootDelta)) {
    sleep(1);
    DivRem(q, r, b0 + FloorRootDelta, 2*a0);

    if(debug) {
      std::cout << "FloorRootDelta is " << FloorRootDelta << std::endl;
      std::cout << "a is " << a0 << std::endl;
      std::cout << "b is " << b0 << std::endl;
      std::cout << "c is " << c0 << std::endl;
      std::cout << "q is " << q << std::endl;
      std::cout << "r is " << r << std::endl;
    }

    c1 = -a0;
    b1 = FloorRootDelta - r;
    a1 = q*((b0 - b1)/2) - c0;

    BTemp = B1;
    B1 = q*B1 + B0;
    B0 = B1;

    a0 = a1;
    b0 = b1;
    c0 = c1;

    if(debug) {
      std::cout << "a1 is " << a1 << std::endl;
      std::cout << "b1 is " << b1 << std::endl;
      std::cout << "c1 is " << c1 << std::endl;

      std::cout << "B0 is " << B0 << std::endl;
      std::cout << "B1 is " << B1 << std::endl;
    }
  }


  RelativeGenerator->set_abd(2*B0*a0 + B1*b0, B1, 2*a0);
  //RelativeGenerator->invert();

  if(RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, ZZ(-1));
  }

  if(debug) {
    std::cout << "RG is " << RelativeGenerator->conv_RR() << std::endl;
    std::cout << "distance is " << log(RelativeGenerator->conv_RR()) << std::endl;
  }


  A.assign(a0, b0, c0);
}
