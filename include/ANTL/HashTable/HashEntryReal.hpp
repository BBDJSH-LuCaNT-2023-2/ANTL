/**
 * @file HashEntryReal.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_QO_HASH_ENTRY_REAL_H
#define ANTL_QUADRATIC_QO_HASH_ENTRY_REAL_H

#include <ANTL/HashTable/HashEntry.hpp>

// Currently using NTL::RR in place of qo_distance; all changes marked by
// comment #include <ANTL/quadratic/qo_distance.hpp>

namespace ANTL {
template <class T> class HashEntryReal;

template <class T>
std::istream &operator>>(std::istream &in, HashEntryReal<T> &A);

template <class T>
std::ostream &operator<<(std::ostream &out, const HashEntryReal<T> &A);

//
// Class: HashEntryReal<T>
//
// This class represents one element in the hash table.  It is essentially
//      an element in a linked list.
//

template <class T> class HashEntryReal : public HashEntry<T> {
protected:
  RR d;
  // qo_distance < T > d;

public:
  HashEntryReal(){};

  HashEntryReal(const T &na, const T &nb, const RR &nd);
  // HashEntryReal(const T &na, const T &nb, const qo_distance<T> &nd);

  ~HashEntryReal(){};

  RR get_d() const { return d; };
  // qo_distance<T> get_d() const { return d; };

  //
  // input/output
  //

  friend std::istream &operator>><T>(std::istream &in, HashEntryReal<T> &A);
  friend std::ostream &operator<<<T>(std::ostream &out,
                                     const HashEntryReal<T> &A);
};

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../src/HashTable/HashEntryReal_impl.hpp"

#endif // guard
