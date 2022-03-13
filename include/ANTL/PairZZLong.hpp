/**
 * @file PairZZLong.hpp
 * @author Mark Velichka
 * @version $Header$
 *
 * A new PairZZLong class, inherrited from a PairZZLong created by NTL's
 * macros. Adds a hash_value() funtion to the class.
 */

#ifndef ANTL_NTL_PairZZLong__H
#define ANTL_NTL_PairZZLong__H

#include <NTL/ZZ.h>
#include <NTL/pair.h>

NTL_CLIENT

namespace ANTL {
//
// PairZZLong_base is an NTL-compatible class.
//
/*
  NTL_pair_decl (ZZ, long, PairZZLong_base)
  NTL_pair_io_decl (ZZ, long, PairZZLong_base)
  NTL_pair_eq_decl (ZZ, long, PairZZLong_base)
*/

/**
 * @brief An extension of an NTL ZZ-long pair with support for hashing.
 */
//    class PairZZLong : public PairZZLong_base
class PairZZLong : public Pair<ZZ, long> {
public:
  // constructors
  // default
  PairZZLong() : Pair<ZZ, long>() {
#ifdef DEBUG_2
    cout << "default: a = " << a << " b = " << b << endl;
#endif
  }

  // pre-initialized
  PairZZLong(const ZZ &a, const long b) {
    this->a = a;
    this->b = b;
#ifdef DEBUG_2
    cout << "overload: a = " << a << " b = " << b << endl;
    cout << this->a << " " << this->b << endl;
#endif
  }

  // copy
  PairZZLong(const PairZZLong &in) {
    a = in.a;
    b = in.b;
#ifdef DEBUG_2
    cout << "Overload2: a = " << a << " b = " << b << endl;
#endif
  }

  // overloading =
  PairZZLong &operator=(const PairZZLong &in) {
    a = in.a;
    b = in.b;
#ifdef DEBUG_2
    cout << "op: in.a = " << in.a << " in.b = " << in.b;
    cout << endl;
    cout << "op: a = " << a << " b = " << b << endl;
#endif
    return *this;
  }

  /// get the ZZ
  ZZ get_a() const { return a; }

  /// get the long
  long get_b() const { return b; }
};

/**
 * @brief Returns it's hash value.
 */
ZZ hash_value(const PairZZLong &A);

} // namespace ANTL

#endif // guard
