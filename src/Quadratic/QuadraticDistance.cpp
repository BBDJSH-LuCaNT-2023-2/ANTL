/**
 * @file QuadraticDistance.cpp
 * @author Jonathan Hammell
 * @version $Header$
 */

#include <ANTL/Quadratic/QuadraticDistance.hpp>

namespace ANTL {

/*
  template < class T >
  void add(QuadraticDistance<T> &c, const QuadraticDistance<T> &a, const
  QuadraticDistance<T> &b)
  {
  NTL::add(c,a,b);
  }

  template < class T >
  void add(QuadraticDistance<T> &c, const QuadraticDistance<T> &a, const NTL::ZZ
  &b)
  {
  NTL::add(c,a,b);
  }

  template < class T >
  void subtract(QuadraticDistance<T> &c, const QuadraticDistance<T> &a, const
  QuadraticDistance<T> &b)
  {
  NTL::sub(c,a,b);
  }

  template < class T >
  void subtract(QuadraticDistance<T> &c, const QuadraticDistance<T> &a, const
  NTL::ZZ &b)
  {
  NTL::sub(c,a,b);
  }
*/

template <>
void add<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
             const QuadraticDistance<ZZ> &b) {
  multiply((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
           (QuadraticDistanceBase &)b);
}

template <>
void add<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
             const NTL::ZZ &b) {
  multiply((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
           (QuadraticDistanceBase &)b);
}

template <>
void subtract<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
                  const QuadraticDistance<ZZ> &b) {
  divide((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
         (QuadraticDistanceBase &)b);
}

template <>
void subtract<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
                  const NTL::ZZ &b) {
  divide((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
         (QuadraticDistanceBase &)b);
}

template <>
void add<long>(QuadraticDistance<long> &c, const QuadraticDistance<long> &a,
               const QuadraticDistance<long> &b) {
  multiply((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
           (QuadraticDistanceBase &)b);
}

template <>
void add<long>(QuadraticDistance<long> &c, const QuadraticDistance<long> &a,
               const NTL::ZZ &b) {
  multiply((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
           (QuadraticDistanceBase &)b);
}

template <>
void subtract<long>(QuadraticDistance<long> &c,
                    const QuadraticDistance<long> &a,
                    const QuadraticDistance<long> &b) {
  divide((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
         (QuadraticDistanceBase &)b);
}

template <>
void subtract<long>(QuadraticDistance<long> &c,
                    const QuadraticDistance<long> &a, const NTL::ZZ &b) {
  divide((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
         (QuadraticDistanceBase &)b);
}

template <>
void add<long long>(QuadraticDistance<long long> &c,
                    const QuadraticDistance<long long> &a,
                    const QuadraticDistance<long long> &b) {
  multiply((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
           (QuadraticDistanceBase &)b);
}

template <>
void add<long long>(QuadraticDistance<long long> &c,
                    const QuadraticDistance<long long> &a, const NTL::ZZ &b) {
  multiply((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
           (QuadraticDistanceBase &)b);
}

template <>
void subtract<long long>(QuadraticDistance<long long> &c,
                         const QuadraticDistance<long long> &a,
                         const QuadraticDistance<long long> &b) {
  divide((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
         (QuadraticDistanceBase &)b);
}

template <>
void subtract<long long>(QuadraticDistance<long long> &c,
                         const QuadraticDistance<long long> &a,
                         const NTL::ZZ &b) {
  divide((QuadraticDistanceBase &)c, (QuadraticDistanceBase &)a,
         (QuadraticDistanceBase &)b);
}

} // namespace ANTL
