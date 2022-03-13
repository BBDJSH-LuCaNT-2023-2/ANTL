/**
 * @file IndexedHashTable.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_INDEXED_HASH_TABLE_H
#define ANTL_INDEXED_HASH_TABLE_H

#include <ANTL/HashTable/HashTable.hpp>

namespace ANTL {

//
// Class: IndexedHashTable<T>
//
// This class represents a hash table, whose elemets are of some type T.  The
// elements can be accessed either as in a regular hash table or sequentially,
// as in a list.
//
template <class T> class IndexedHashTable : public HashTable<T> {
  using HashTable<T>::size;
  using HashTable<T>::curr_size;
  using HashTable<T>::buckets;
  using HashTable<T>::last_one;

private:
  long allocated; // size of array

  HEntry<T> **IDX; // array of pointers to list elements
  long outstyle;   // output style (0=list, 1=hash table)

public:
  IndexedHashTable();
  ~IndexedHashTable();

  void assign(const IndexedHashTable<T> &old_table);
  IndexedHashTable<T> &operator=(const IndexedHashTable<T> &old_table);

  const T operator[](long i) const;
  const T member(long i) const;
  T &member_ref(long i);

  void remove(const T &G);
  void remove_from(long i);
  void empty();
  void hash(const T &G);

  void output_style(int style);
  void read(std::istream &in);
  void print(std::ostream &out) const;

  friend std::istream &operator>>(std::istream &in, IndexedHashTable<T> &HT) {
    HT.read(in);
    return (in);
  }

  friend std::ostream &operator<<(std::ostream &out,
                                  const IndexedHashTable<T> &HT) {
    HT.print(out);
    return (out);
  }
};

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../src/HashTable/IndexedHashTable_impl.hpp"

#endif
