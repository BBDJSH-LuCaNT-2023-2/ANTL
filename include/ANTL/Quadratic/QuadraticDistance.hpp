/**
 * @file QuadraticDistance.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_QO_DISTANCE_H
#define ANTL_QUADRATIC_QO_DISTANCE_H

#include <ANTL/Quadratic/QuadraticDistanceBase.hpp>
#include <NTL/ZZ.h>

namespace ANTL {

template <class T> class QuadraticDistance;

template <class T>
void add(QuadraticDistance<T> &c, const QuadraticDistance<T> &a,
         const QuadraticDistance<T> &b) {
  NTL::add(c, a, b);
};
template <class T>
void add(QuadraticDistance<T> &c, const QuadraticDistance<T> &a,
         const NTL::ZZ &b) {
  NTL::add(c, a, b);
};
template <class T>
void subtract(QuadraticDistance<T> &c, const QuadraticDistance<T> &a,
              const QuadraticDistance<T> &b) {
  NTL::sub(c, a, b);
};
template <class T>
void subtract(QuadraticDistance<T> &c, const QuadraticDistance<T> &a,
              const NTL::ZZ &b) {
  NTL::sub(c, a, b);
};

//
// Class: QuadraticDistance
//
// This class represents a distance (used in infrastructure computations).  For
// number fields, this is a QuadraticDistanceBase.  For function fields, this is
// a ZZ.
//
template <class T> class QuadraticDistance : public ZZ {
public:
  QuadraticDistance() { clear(*this); };
  QuadraticDistance(const ZZ &a) : ZZ(a){};
  void init_distance() { clear(*this); };
  void assign_one() { set(*this); }
  void assign(const QuadraticDistance<T> &dist) { *this = dist; }
  ZZ eval() const { return *this; };

  ZZ get_elog() const { return *this; };

  void adjust_inverse(const T &a) { *this -= deg(a); };

  friend void add<T>(QuadraticDistance<T> &c, const QuadraticDistance<T> &a,
                     const QuadraticDistance<T> &b);
  friend void add<T>(QuadraticDistance<T> &c, const QuadraticDistance<T> &a,
                     const NTL::ZZ &b);
  friend void subtract<T>(QuadraticDistance<T> &c,
                          const QuadraticDistance<T> &a,
                          const QuadraticDistance<T> &b);
  friend void subtract<T>(QuadraticDistance<T> &c,
                          const QuadraticDistance<T> &a, const NTL::ZZ &b);
};

template <> class QuadraticDistance<ZZ> : public QuadraticDistanceBase {
public:
  QuadraticDistance() : QuadraticDistanceBase() {}
  QuadraticDistance(const ZZ &nd, const ZZ &nk)
      : QuadraticDistanceBase(nd, nk) {}

  void init_distance() { assign_one(); };
  void assign(const QuadraticDistance<ZZ> &dist) {
    d = dist.get_d();
    k = dist.get_k();
  };
  ZZ eval() const { return k; };
  ZZ get_elog() const { return CeilToZZ(get_log()); };

  void adjust_inverse(const ZZ &a) { this->divide(a); };

  friend void add<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
                      const QuadraticDistance<ZZ> &b);
  friend void add<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
                      const NTL::ZZ &b);
  friend void subtract<ZZ>(QuadraticDistance<ZZ> &c,
                           const QuadraticDistance<ZZ> &a,
                           const QuadraticDistance<ZZ> &b);
  friend void subtract<ZZ>(QuadraticDistance<ZZ> &c,
                           const QuadraticDistance<ZZ> &a, const NTL::ZZ &b);
};

template <> class QuadraticDistance<long> : public QuadraticDistanceBase {
public:
  void init_distance() { assign_one(); };
  void assign(const QuadraticDistance<long> &dist) {
    d = dist.get_d();
    k = dist.get_k();
  };
  ZZ eval() const { return k; };

  ZZ get_elog() const { return CeilToZZ(get_log()); };

  void adjust_inverse(const long &a) { this->divide(a); };

  friend void add<long>(QuadraticDistance<long> &c,
                        const QuadraticDistance<long> &a,
                        const QuadraticDistance<long> &b);
  friend void add<long>(QuadraticDistance<long> &c,
                        const QuadraticDistance<long> &a, const NTL::ZZ &b);
  friend void subtract<long>(QuadraticDistance<long> &c,
                             const QuadraticDistance<long> &a,
                             const QuadraticDistance<long> &b);
  friend void subtract<long>(QuadraticDistance<long> &c,
                             const QuadraticDistance<long> &a,
                             const NTL::ZZ &b);
};

template <> class QuadraticDistance<long long> : public QuadraticDistanceBase {
public:
  void init_distance() { assign_one(); };
  void assign(const QuadraticDistance<long long> &dist) {
    d = dist.get_d();
    k = dist.get_k();
  };
  ZZ eval() const { return k; };

  ZZ get_elog() const { return CeilToZZ(get_log()); };

  void adjust_inverse(const long long &a) { this->divide(a); };

  friend void add<long long>(QuadraticDistance<long long> &c,
                             const QuadraticDistance<long long> &a,
                             const QuadraticDistance<long long> &b);
  friend void add<long long>(QuadraticDistance<long long> &c,
                             const QuadraticDistance<long long> &a,
                             const NTL::ZZ &b);
  friend void subtract<long long>(QuadraticDistance<long long> &c,
                                  const QuadraticDistance<long long> &a,
                                  const QuadraticDistance<long long> &b);
  friend void subtract<long long>(QuadraticDistance<long long> &c,
                                  const QuadraticDistance<long long> &a,
                                  const NTL::ZZ &b);
};

template <>
void add<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
             const QuadraticDistance<ZZ> &b);
template <>
void add<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
             const NTL::ZZ &b);
template <>
void subtract<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
                  const QuadraticDistance<ZZ> &b);
template <>
void subtract<ZZ>(QuadraticDistance<ZZ> &c, const QuadraticDistance<ZZ> &a,
                  const NTL::ZZ &b);
template <>
void add<long>(QuadraticDistance<long> &c, const QuadraticDistance<long> &a,
               const QuadraticDistance<long> &b);
template <>
void add<long>(QuadraticDistance<long> &c, const QuadraticDistance<long> &a,
               const NTL::ZZ &b);
template <>
void subtract<long>(QuadraticDistance<long> &c,
                    const QuadraticDistance<long> &a,
                    const QuadraticDistance<long> &b);
template <>
void subtract<long>(QuadraticDistance<long> &c,
                    const QuadraticDistance<long> &a, const NTL::ZZ &b);
template <>
void add<long long>(QuadraticDistance<long long> &c,
                    const QuadraticDistance<long long> &a,
                    const QuadraticDistance<long long> &b);
template <>
void add<long long>(QuadraticDistance<long long> &c,
                    const QuadraticDistance<long long> &a, const NTL::ZZ &b);
template <>
void subtract<long long>(QuadraticDistance<long long> &c,
                         const QuadraticDistance<long long> &a,
                         const QuadraticDistance<long long> &b);
template <>
void subtract<long long>(QuadraticDistance<long long> &c,
                         const QuadraticDistance<long long> &a,
                         const NTL::ZZ &b);

} // namespace ANTL

#endif // guard
