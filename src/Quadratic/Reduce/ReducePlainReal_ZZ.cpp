/**
 * @file qo_reduce_plain_ZZ.cpp
 * @author Michael Jacobson
 * @remark Basic ideal reduction specialization (ZZ base type).
 */

#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>

// ReducePlainReal<ZZ>::Reduce
// Task: Reduces the ideal

template <> void ReducePlainReal<ZZ>::reduce(QuadraticIdealBase<ZZ> &A) {

  ZZ a_old = A.get_a(), b_old = A.get_b(), c_old = A.get_c();
  ZZ a_new, b_new, c_new, q, r, B_m1 = ZZ(0), B_m2 = ZZ(1), B_temp;

  set(*RelativeGenerator);

  if (!A.is_normal()) {
    A.normalize();
  }

  if (A.is_reduced()) {
    return;
  }

  do {
    DivRem(q, r, b_old + FloorRootDelta, 2 * a_old);

    c_new = -a_old;
    b_new = FloorRootDelta - r;
    a_new = q * ((b_old - b_new) / 2) - c_old;

    B_temp = B_m1;
    B_m1 = q * B_m1 + B_m2;
    B_m2 = B_temp;

    a_old = a_new;
    b_old = b_new;
    c_old = c_new;

    abs(a_new, a_new);

    if(FloorRootDelta < a_new) {
      b_new += ((a_new-b_new) / (2*a_new)) * (2*a_new);
    }

    else {
      b_new += ((FloorRootDelta-b_new) / (2*a_new)) * (2*a_new);
    }

    c_new = (Delta - (b_new*b_new))/(-4*a_new);

  } while (!((abs(FloorRootDelta - 2 * a_new) < b_new) && (b_new < FloorRootDelta + 1)));

  RelativeGenerator->set_abd(2 * B_m2 * a_old + B_m1 * b_old, -B_m1, 2 * a_old);
  RelativeGenerator->invert();

  if (RelativeGenerator->conv_RR() < 0) {
    mul(*RelativeGenerator, *RelativeGenerator, ZZ(-1));
  }

  A.assign(a_new, b_new, c_new);
}
