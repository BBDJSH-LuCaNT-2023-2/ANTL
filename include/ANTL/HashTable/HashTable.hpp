/**
 * @file HashTable.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_HASH_TABLE_H
#define ANTL_HASH_TABLE_H

#include <ANTL/common.hpp>

namespace ANTL {

/**
 * This class represents one element in the hash table.  It is essentially
 * an element in a linked list.
 */
template <class T> class HEntry {
protected:
public:
  T *item;         // pointer to the data item
  HEntry<T> *next; // pointers to next and previous list elements

  HEntry() {
    item = NULL;
    next = NULL;
  };
  ~HEntry();
};

//
// Class: HashTable<T>
//
// This class represents a hash table, whose elements are of some type T.
//	The collision resolution scheme is bucketing, and each bucket is a
//	linked list.
// The hash function is simply the key value of the data item modulo the
//	number of buckets in the table.  The user must define a key function
//	corresponding to the type T, and initialize it with the function
//	set_key_function().
//

template <class T> class HashTable {
protected:
  long size;           // number of buckets in the hash table
  long curr_size;      // current number of elemets
  HEntry<T> **buckets; // array of data buckets
  T *last_one;         // pointer to most recent entry

public:
  HashTable();
  ~HashTable();

  void assign(const HashTable<T> &old_table);
  HashTable<T> &operator=(const HashTable<T> &old_table);

  void initialize(const long table_size);

  long no_of_buckets() const;
  long no_of_elements() const;

  void remove(const T &G);
  void empty();
  const T last_entry() const;
  void hash(const T &G);
  T *search(const T &G)
#if (defined(__GNUC__) || !defined(__mips) || !defined(__sgi) || defined(__EDG))
      const
#endif
      ;

  HEntry<T> *get_bucket(const T &G)
#if (defined(__GNUC__) || !defined(__mips) || !defined(__sgi) || defined(__EDG))
      const
#endif
      ;

  void read(std::istream &in);
  void print(std::ostream &out) const;
  void write(std::ostream &out) const;

  friend std::istream &operator>>(std::istream &in, HashTable<T> &HT) {
    HT.read(in);
    return (in);
  }

  friend std::ostream &operator<<(std::ostream &out, const HashTable<T> &HT) {
    HT.print(out);
    return (out);
  }
};

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../src/HashTable/HashTable_impl.hpp"

#endif // guard
