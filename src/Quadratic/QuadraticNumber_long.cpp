/**
 * @file QuadraticNumber_long.cpp
 * @author Michael Jacobson
 * @remark Quadratic number function specializations (long base type).
 */

#include <ANTL/Quadratic/QuadraticNumber.hpp>

template <> void QuadraticNumber<long>::normalize() {
  if (d < 0) {
    a = -a;
    b = -b;
    d = -d;
  }
  long g = GCD(GCD(a, b), d);
  if (g != 1) {
    a /= g;
    b /= g;
    d /= g;
  }
}

template <> bool QuadraticNumber<long>::isUnit() const {
  QQ<long> tempQ = QQ<long>();
  abs(tempQ, getNorm());
  return tempQ == 1;
}
