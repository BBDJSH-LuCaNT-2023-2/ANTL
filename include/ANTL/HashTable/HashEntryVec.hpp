/**
 * @file HashEntryVec.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_HASH_ENTRY_VEC_H
#define ANTL_QUADRATIC_HASH_ENTRY_VEC_H

#include <ANTL/HashTable/HashEntry.hpp>
#include <NTL/vec_ZZ.h>

namespace ANTL {
template <class T> class HashEntryVec;

template <class T>
std::istream &operator>>(std::istream &in, HashEntryVec<T> &A);

template <class T>
std::ostream &operator<<(std::ostream &out, const HashEntryVec<T> &A);

//
// Class: qo_hash_hentry_vec<T>
//
// This class represents one element in the hash table.  It is essentially
//      an element in a linked list.
//

template <class T> class HashEntryVec : public HashEntry<T> {
protected:
  vec_ZZ d;

public:
  HashEntryVec(){};
  HashEntryVec(const T &na, const T &nb, const vec_ZZ &newd);

  ~HashEntryVec(){};

  vec_ZZ get_d() const { return d; };

  //
  // input/output
  //

  friend std::istream &operator>><T>(std::istream &in, HashEntryVec<T> &A);
  friend std::ostream &operator<<<T>(std::ostream &out,
                                     const HashEntryVec<T> &A);
};

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../src/HashTable/HashEntryVec_impl.hpp"

#endif // guard
