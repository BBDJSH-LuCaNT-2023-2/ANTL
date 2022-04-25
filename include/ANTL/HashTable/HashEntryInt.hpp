/**
 * @file HashEntryInt.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_HASH_ENTRY_INT_H
#define ANTL_QUADRATIC_HASH_ENTRY_INT_H

#include <ANTL/HashTable/HashEntry.hpp>
#include <NTL/ZZ.h>

namespace ANTL {

template <class T, class S> class HashEntryInt;

template <class T, class S>
std::istream &operator>>(std::istream &in, HashEntryInt<T, S> &A);

template <class T, class S>
std::ostream &operator<<(std::ostream &out, const HashEntryInt<T, S> &A);

//
// Class: qo_hash_hentryInt<T,S>
//
// This class represents one element in the hash table.  It is essentially
//      an element in a linked list.
//

template <class T, class S> class HashEntryInt : public HashEntry<T> {
protected:
  S d;

public:
  HashEntryInt(){};
  HashEntryInt(const T &na, const T &nb, const S &newd);

  ~HashEntryInt(){};

  S get_d() const { return d; };

  //
  // input/output
  //

  friend std::istream &operator>>
      <T, S>(std::istream &in, HashEntryInt<T, S> &A);
  friend std::ostream &operator<<<T, S>(std::ostream &out,
                                        const HashEntryInt<T, S> &A);
};

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../src/HashTable/HashEntryInt_impl.hpp"

#endif // guard
