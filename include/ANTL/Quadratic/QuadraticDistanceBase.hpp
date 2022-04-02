/**
 * @file QuadraticDistanceBase.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_QO_DISTANCE_BASE_H
#define ANTL_QUADRATIC_QO_DISTANCE_BASE_H

//#include <strstream>
//#include <strstream.h>

#include <ANTL/common.hpp>
#include <NTL/RR.h>
#include <NTL/ZZ.h>
#include <stdio.h>
#include <string>

namespace ANTL {
//
// Class: QuadraticDistanceBase
//
// This class represents an extended-fp-approximation of a principal ideal
// generator. It is used for keeping track of distances of ideals in real
// quadratic fields (i.e., for infrastructure computations).
//
// It consists of a pair (d, k) such that 2^k d/2^p approximates the generator.
// p is a global precision constant.
//
// Before using any QuadraticDistanceBase objects, the static function
// QuadraticDistanceBase::set_precision must be called to set the precision
// constants and to compute approximations of the square root of the current
// quadratic field discriminant.  This is done automatically whenever the static
// set_current_order function of the qi_pair class is called.
//
// Two output modes are possible.  In QuadraticDistanceBase::INTEGER mode, the
// QuadraticDistanceBase is output as an ordered pair (d, k).  In
// QuadraticDistanceBase::REAL mode (the default), an approximation of the
// natural logarithm of the generator is computed and output.  To use the input
// operator, INTEGER mode must be set.
//
#define QD_REAL 0
#define QD_INTEGER 1
class QuadraticDistanceBase {
protected:
  // approximation of the principal ideal generator
  ZZ d;
  ZZ k;

  // approximations of sqrt(Delta) to various precision
  static ZZ rootDs;
  static ZZ rootDp;
  static ZZ rootDpm;

  // precision constants
  static long s;
  static long m;
  static long p;

  // powers of 2
  static ZZ s2;
  static ZZ m2;
  static ZZ p2;  // 2^p
  static ZZ p21; // 2^(p + 1)
  static ZZ p22; // 2^(2p)
  static RR log2;

  static ZZ xi;
  static ZZ yi;

  // output mode
  static long output_mode;

public:
  static const long REAL = 0;
  static const long INTEGER = 1;

  QuadraticDistanceBase(){};
  QuadraticDistanceBase(const ZZ &nd, const ZZ &nk);
  QuadraticDistanceBase(const QuadraticDistanceBase &qd);

  ~QuadraticDistanceBase(){};

  static void set_precision(const ZZ &Delta);
  static void set_precision(const ZZ &Delta, const ZZ &S);
  static void set_precision(long Delta);
  static void set_precision(long Delta, const ZZ &S);
  static void set_precision(long long Delta);
  static void set_precision(long long Delta, const ZZ &S);
  static void set_output_mode(long mode);

  static long get_s() { return s; };
  static long get_m() { return m; };
  static long get_p() {
    //	  cout << "Doesn't work!" << endl;
    //	  return 0;
    return p;
  };
  static long get_output_mode() { return output_mode; };

  static ZZ get_rootDs() { return rootDs; };
  static ZZ get_rootDp() { return rootDp; };
  static ZZ get_rootDpm() { return rootDpm; };

  static ZZ get_xi() { return xi; };
  static ZZ get_yi() { return yi; };

  void assign_zero();
  void assign_one();
  void assign(const ZZ &d2, const ZZ &k2) {
    d = d2;
    k = k2;
  };
  void assign(const QuadraticDistanceBase &qd);
  QuadraticDistanceBase &operator=(const QuadraticDistanceBase &qd);

  friend void clear(QuadraticDistanceBase &qd);
  friend void set(QuadraticDistanceBase &qd);

  ZZ get_d() const { return d; };
  ZZ get_k() const { return k; };
  RR get_log() const {
    if (is_zero())
      return to_RR(0);
    else {
      RR logqd = to_RR(k - p) * log2 + log(to_RR(d));
      return logqd;
    }
  };

  bool is_zero() const { return ((d == p2 || IsZero(d)) && IsZero(k)); };
  bool is_one() const { return IsZero(k) && d == p2; };

  friend bool IsZero(const QuadraticDistanceBase &qd) { return qd.is_zero(); };

  friend bool operator==(const QuadraticDistanceBase &A,
                         const QuadraticDistanceBase &B);
  friend bool operator!=(const QuadraticDistanceBase &A,
                         const QuadraticDistanceBase &B);

  friend void swap(QuadraticDistanceBase &A, QuadraticDistanceBase &B);

  void normalize();

  void invert();

  void multiply(const ZZ &x);
  void multiply(long x);
  void multiply(long long x);

  void multiply_reduce(const ZZ &x, const ZZ &y);
  void multiply_reduce(long x, long y);
  void multiply_reduce(long long x, long long y);

  void multiply_rho(const ZZ &x, const ZZ &y);
  void multiply_rho(long x, long y);
  void multiply_rho(long long x, long long y);

  void multiply_inverse_rho(const ZZ &q, const ZZ &x, const ZZ &y);
  void multiply_inverse_rho(long q, long x, long y);
  void multiply_inverse_rho(long long q, long long x, long long y);

  void divide(const ZZ &x);
  void divide(long x);
  void divide(long long x);

  void reduce_div(const ZZ &x);

  friend void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
                       const ZZ &x, const ZZ &y);
  friend void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
                       const long &x, const long &y);
  friend void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
                       const long long &x, const long long &y);

  friend void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
                       const QuadraticDistanceBase &B);
  friend void square(QuadraticDistanceBase &C, const QuadraticDistanceBase &A);
  friend void divide(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
                     const QuadraticDistanceBase &B);

  friend std::istream &operator>>(std::istream &in, QuadraticDistanceBase &A);
  friend std::ostream &operator<<(std::ostream &out,
                                  const QuadraticDistanceBase &A);
};

} // namespace ANTL

#endif // guard
