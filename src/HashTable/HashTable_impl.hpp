/**
 * @file HashTable_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from HashTable.h only.
 * @version $Header$
 */

#include <string.h>

namespace ANTL {

//
// HEntry destructor:
//

template <class T> HEntry<T>::~HEntry() {
  if (item)
    delete item;
}

//
// constructor:
//	- set number of buckets and number of elements to 0
//	- set all polongers to NULL
//

template <class T> HashTable<T>::HashTable() {
  size = 0;
  curr_size = 0;
  buckets = (HEntry<T> **)NULL;
  last_one = (T *)NULL;
}

//
// destructor:
//	if the table has been allocated, delete each bucket with the function
//	empty(), and then delete the list of buckets
//

template <class T> HashTable<T>::~HashTable() {
  if (buckets) {
    empty();
    delete[] buckets;
    buckets = (HEntry<T> **)NULL;
    last_one = (T *)NULL;
  }
}

//
// HashTable<T>::assign()
//
// Task:
//	make a copy of an existing hash table
//

template <class T> void HashTable<T>::assign(const HashTable<T> &old_table) {
  long i;
  HEntry<T> *ptr, *newone;
  T *tptr;

  // initialize new table
  initialize(old_table.size);

  // perform linked list traversal, and hash each element longo the new table
  newone = NULL;
  for (i = 0; i < old_table.size; ++i) {
    ptr = old_table.buckets[i];
    while (ptr) {
      tptr = ptr->item;
      hash(*tptr);
      ptr = ptr->next;
    }
  }
}

//
// operator =
//
// Task:
//	make a copy of an existing hash table
//

template <class T>
HashTable<T> &HashTable<T>::operator=(const HashTable<T> &old_table) {
  this->assign(old_table);
  return *this;
}

//
// HashTable<T>::initialize()
//
// Task:
//	set the number of buckets and allocate storage for the array of buckets
////

template <class T> void HashTable<T>::initialize(const long table_size) {
  ZZ temp;

  if (buckets) {
    empty();
    delete[] buckets;
    buckets = (HEntry<T> **)NULL;
    last_one = (T *)NULL;
  }

  // number of buckets should be a prime
  temp = NextPrime(table_size - 1);
  conv(size, temp);

  buckets = new HEntry<T> *[size];
  memset(buckets, '\0', (long)(size) * (sizeof buckets[0]));
}

//
// HashTable<T>::no_of_buckets()
//
// Task:
//	returns the number of buckets in the table
//

template <class T> long HashTable<T>::no_of_buckets() const { return size; }

//
// HashTable<T>::no_of_elements()
//
// Task:
//	returns the number of elements currently in the table
//

template <class T> long HashTable<T>::no_of_elements() const {
  return curr_size;
}

//
// HashTable<T>::remove()
//
// Task:
//	removes G from the table, if it is there.
//

template <class T> void HashTable<T>::remove(const T &G) {
  T *target, *tptr;
  HEntry<T> *ptr, *pptr, *nptr;
  long i, j;

  // check whether G is in the table (linked list search)
  target = (T *)NULL;
  i = hash_value(G) % size;
  if (i < 0)
    i += size;
  ptr = buckets[i];
  pptr = NULL;
  while ((ptr) && (!target)) {
    tptr = ptr->item;
    if (G == *tptr)
      target = tptr;
    else {
      pptr = ptr;
      ptr = ptr->next;
    }
  }

  // if so, remove it (linked list delete)
  if (target) {
    nptr = ptr->next;

    if (pptr)
      pptr->next = nptr;
    else
      // G was the first element in the list - handle specially
      buckets[i] = nptr;

    delete ptr;
    --curr_size;
  }
}

//
// HashTable<T>::empty()
//
// Task:
//	delete all the elements in the table.  At the end, the number of
//	buckets is the same, but they will all be empty.
//

template <class T> void HashTable<T>::empty() {
  long i;
  HEntry<T> *ptr, *nptr;

  // execute a linked list delete on each bucket
  for (i = 0; i < size; ++i) {
    ptr = buckets[i];
    while (ptr) {
      nptr = ptr->next;
      delete ptr;
      ptr = nptr;
    }
    buckets[i] = (HEntry<T> *)NULL;
  }

  curr_size = 0;
  last_one = (T *)NULL;
}

//
// HashTable<T>::last_entry()
//
// Task:
//	returns a constant reference to the most recent element inserted in
//	the hash table
//
// Conditions:
//	there must be at least one element in the list
//

template <class T> const T HashTable<T>::last_entry() const {
  if (!last_one) {
    cerr << "HashTable::last_entry - table is empty" << endl;
    exit(1);
  }

  return *last_one;
}

//
// HashTable<T>::hash()
//
// Task:
//	insert G into the hash table
//

template <class T> void HashTable<T>::hash(const T &G) {
  long i;
  HEntry<T> *ptr, *nptr, *newone;
  T *newT;

  // allocate new list element and a copy of G
  newone = new HEntry<T>;
  newT = new T;

  (*newT) = G;
  newone->item = newT;

  ++curr_size;
  last_one = newT;

  // compute correct bucket and append G to the end of the linked list
  i = hash_value(G) % size;
  if (i < 0)
    i += size;

  ptr = buckets[i];
  if (ptr) {
    while (ptr) {
      nptr = ptr;
      ptr = nptr->next;
    }
    nptr->next = newone;
  } else
    buckets[i] = newone;
}

//
// HashTable<T>::search()
//
// Task:
//	returns a polonger to the first occurance of G in the hash table.  If
//	G is not in the hash table, the NULL polonger is returned.
//

template <class T>
T *HashTable<T>::search(const T &G)
#if (defined(__GNUC__) || !defined(__mips) || !defined(__sgi) || defined(__EDG))
    const
#endif

{
  long i;
  T *target, *tptr;
  HEntry<T> *ptr;

  // compute correct bucket and perform linked list search
  i = hash_value(G) % size;
  if (i < 0)
    i += size;
  target = (T *)NULL;
  ptr = buckets[i];
  while ((ptr) && (!target)) {
    tptr = ptr->item;
    //    cout << "tptr = " << tptr << endl;
    //    cout << "(*tptr).number = " << (*tptr).number << endl;
    if (G == (*tptr)) {
//       std::cout << "found G == *tptr!" << std::endl;
      target = tptr;
    }
    else
      ptr = ptr->next;
  }

  return target;
}

//
// HashTable<T>::get_bucket()
//
// Task:
//      returns a polonger to the bucket corresponding to the element G.
//      This polonger is essentially the head of a simply linked list.
//      Note that this function does not test if the bucket is full or
//      empty - this is the user's responsibility.
//

template <class T>
HEntry<T> *HashTable<T>::get_bucket(const T &G)
#if (defined(__GNUC__) || !defined(__mips) || !defined(__sgi) || defined(__EDG))
    const
#endif
{
  long i;

  // compute correct bucket index
  i = hash_value(G) % size;
  if (i < 0)
    i += size;

  return buckets[i];
}

//
// HashTable<T>::read()
//
// Task:
//	read in a hash table from the istream in.
//
// Conditions:
//	input must be the number of buckets on one line, followed by the
//	number of elements to insert longo the list on a new line, and finally
//	n instances of type T.
//

template <class T> void HashTable<T>::read(std::istream &in) {
  long i;
  long new_size, num;
  T new_item;

  in >> new_size;
  initialize(new_size);

  in >> num;
  for (i = 0; i < num; ++i) {
    in >> new_item;
    hash(new_item);
  }
}

//
// HashTable<T>::print()
//
// Task:
//	output a hash table to the ostream out.
//

template <class T> void HashTable<T>::print(std::ostream &out) const {
  long i;
  HEntry<T> *ptr;
  T *tptr;

  for (i = 0; i < size; ++i) {
    if (buckets[i]) {
      out << "[" << i << "]";
      ptr = buckets[i];
      while (ptr) {
        tptr = ptr->item;
        out << " : " << (*tptr);
        ptr = ptr->next;
      }
      out << "\n" << std::endl;
    }
  }
}

//
// HashTable<T>::write()
//
// Task:
//	output a hash table to the ostream out in such a way that ::read() will
//	be able to use the output.
//

template <class T> void HashTable<T>::write(std::ostream &out) const {
  long i;
  HEntry<T> *ptr;
  T *tptr;
  out << this->no_of_buckets() << endl;
  out << this->no_of_elements() << endl;
  for (i = 0; i < size; ++i) {
    if (buckets[i]) {
      // out << "[" << i << "]";
      ptr = buckets[i];
      while (ptr) {
        tptr = ptr->item;
        // out << " : " << (*tptr);
        out << (*tptr) << endl;
        ptr = ptr->next;
      }
      //      out << "\n";
    }
  }
}

} // namespace ANTL
