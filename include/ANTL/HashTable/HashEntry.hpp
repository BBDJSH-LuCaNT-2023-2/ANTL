/**
 * @file HashEntry.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_QO_HASH_ENTRY_H
#define ANTL_QUADRATIC_QO_HASH_ENTRY_H

#include <ANTL/common.hpp>

namespace ANTL {
template <class T> class HashEntry;

template <class T>
bool operator==(const HashEntry<T> &A, const HashEntry<T> &B);

template <class T>
bool operator!=(const HashEntry<T> &A, const HashEntry<T> &B);

template <class T> istream &operator>>(istream &in, HashEntry<T> &A);

template <class T> ostream &operator<<(ostream &out, const HashEntry<T> &A);

//
// Class: HashEntry<T>
//
// This class represents one element in the hash table.  It is essentially
//      an element in a linked list.
//

template <class T> class HashEntry {
protected:
  ZZ a;
  ZZ b;

public:
  HashEntry() {}
  HashEntry(const T &na, const T &nb) {
    a = eval_poly(na, CARDINALITY<T>());
    b = eval_poly(nb % na, CARDINALITY<T>());
  };

  ~HashEntry() {}

  T get_a() const {
    T ap;
    get_poly_modq(ap, a, CARDINALITY<T>());
    return ap;
  }
  T get_b() const {
    T bp;
    get_poly_modq(bp, b, CARDINALITY<T>());
    return bp;
  }
  ZZ hval() const { return a; }

  friend bool operator==<T>(const HashEntry<T> &A, const HashEntry<T> &B);

  friend bool operator!=<T>(const HashEntry<T> &A, const HashEntry<T> &B);

  //
  // input/output
  //

  friend std::istream &operator>><T>(std::istream &in, HashEntry<T> &A);
  friend std::ostream &operator<<<T>(std::ostream &out, const HashEntry<T> &A);
};

template <> class HashEntry<ZZ> {
protected:
  ZZ a;
  ZZ b;

public:
  HashEntry() {}
  HashEntry(const ZZ &na, const ZZ &nb) {
    a = na;
    b = nb;
  };

  ~HashEntry() {}

  ZZ get_a() const { return a; }
  ZZ get_b() const { return b; }
  ZZ hval() const { return a; }

  friend bool operator==<ZZ>(const HashEntry<ZZ> &A, const HashEntry<ZZ> &B);

  friend bool operator!=<ZZ>(const HashEntry<ZZ> &A, const HashEntry<ZZ> &B);

  //
  // input/output
  //

  friend std::istream &operator>><ZZ>(std::istream &in, HashEntry<ZZ> &A);
  friend std::ostream &operator<<<ZZ>(std::ostream &out,
                                      const HashEntry<ZZ> &A);
};

template <> class HashEntry<long> {
protected:
  long a;
  long b;

public:
  HashEntry() {}
  HashEntry(const long &na, const long &nb) {
    a = na;
    b = nb;
  };

  ~HashEntry() {}

  long get_a() const { return a; }
  long get_b() const { return b; }
  ZZ hval() const { return to<ZZ>(a); }

  friend bool operator==
      <long>(const HashEntry<long> &A, const HashEntry<long> &B);

  friend bool operator!=
      <long>(const HashEntry<long> &A, const HashEntry<long> &B);

  //
  // input/output
  //

  friend std::istream &operator>><long>(std::istream &in, HashEntry<long> &A);
  friend std::ostream &operator<<<long>(std::ostream &out,
                                        const HashEntry<long> &A);
};

template <> class HashEntry<long long> {
protected:
  long long a;
  long long b;

public:
  HashEntry() {}
  HashEntry(const long long &na, const long long &nb) {
    a = na;
    b = nb;
  };

  ~HashEntry() {}

  long long get_a() const { return a; }
  long long get_b() const { return b; }
  ZZ hval() const { return ZZ(a); }

  friend bool operator==
      <long long>(const HashEntry<long long> &A, const HashEntry<long long> &B);

  friend bool operator!=
      <long long>(const HashEntry<long long> &A, const HashEntry<long long> &B);

  //
  // input/output
  //

  friend std::istream &operator>>
      <long long>(std::istream &in, HashEntry<long long> &A);
  friend std::ostream &operator<<<long long>(std::ostream &out,
                                             const HashEntry<long long> &A);
};

template <> bool operator==<ZZ>(const HashEntry<ZZ> &A, const HashEntry<ZZ> &B);
template <>
bool operator==
    <long long>(const HashEntry<long long> &A, const HashEntry<long long> &B);
template <>
bool operator==<long>(const HashEntry<long> &A, const HashEntry<long> &B);

template <class T>
ZZ hash_value(const HashEntry<T> &A)
//  { return hash_value_base(A.get_a()); }
{
  return A.hval();
}

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../src/HashTable/HashEntry_impl.hpp"

#endif // guard
