/**
 * @file QuadraticNumber_ZZ.cpp
 * @author Michael Jacobson
 * @remark Quadratic number function specializations (ZZ base type).
 */

#include <ANTL/Quadratic/QuadraticNumber.hpp>

template <> void QuadraticNumber<ZZ>::normalize() {
  if (d < 0) {
    NTL::negate(a, a);
    NTL::negate(b, b);
    NTL::negate(d, d);
  }
  ZZ g = GCD(GCD(a, b), d);
  if (g != 1) {
    div(a, a, g);
    div(b, b, g);
    div(d, d, g);
  }
}

template <> bool QuadraticNumber<ZZ>::isUnit() const {
  QQ<ZZ> tempQ = QQ<ZZ>();
  abs(tempQ, getNorm());
  return tempQ == ZZ(1);
}
