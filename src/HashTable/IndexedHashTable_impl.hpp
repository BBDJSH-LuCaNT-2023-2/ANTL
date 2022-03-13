/**
 * @file IndexedHashTable_impl.hpp
 * @file Michael Jacobson
 * @remarks This file is to be included from index_hashtable.h only.
 * @version $Header$
 */

namespace ANTL {

//
// constructor:
//      - execute HashTable constructor
//      - initialize list of pointers (expand mode, 0 elements)
//      - set default output style (hash table mode)
//

template <class T> IndexedHashTable<T>::IndexedHashTable() : HashTable<T>() {
  IDX = new HEntry<T> *[1];
  allocated = 1;
  outstyle = 0;
}

//
// destructor:
//      reset the list of pointers
//

template <class T> IndexedHashTable<T>::~IndexedHashTable() {
  if (IDX)
    delete[] IDX;
}

//
// IndexedHashTable<T>::assign()
//
// Task:
//    make a copy of an existing hash table
//

template <class T>
void IndexedHashTable<T>::assign(const IndexedHashTable<T> &old_table) {
  long i;

  // initialize new table
  this->initialize(old_table.size);

  // insert every element of the old table into the new table
  for (i = 0; i < old_table.curr_size; ++i)
    hash(old_table[i]);
}

//
// operator =
//
// Task:
//      make a copy of an existing hash table
//

template <class T>
IndexedHashTable<T> &
IndexedHashTable<T>::operator=(const IndexedHashTable<T> &old_table) {
  this->assign(old_table);
  return *this;
}

//
// operator []
//
// Task:
//      return a read-only reference to the ith element (sequentially)
//
// Conditions:
//      the index i must satisfy 0 <= i < curr_size
//

template <class T> const T IndexedHashTable<T>::operator[](long i) const {
  HEntry<T> *ptr;
  T *tptr;

  if ((i < 0) || (i >= this->curr_size)) {
    cerr << "IndexedHashTable::operator [] - i >= current size or < 0" << endl;
    exit(1);
  }

  ptr = IDX[i];

  tptr = ptr->item;

  return *tptr;
}

//
// IndexedHashTable<T>::member()
//
// Task:
//      return a read-only reference to the ith element (sequentially)
//
// Conditions:
//      the index i must satisfy 0 <= i < curr_size
//

template <class T> const T IndexedHashTable<T>::member(long i) const {
  HEntry<T> *ptr;
  T *tptr;

  if ((i < 0) || (i >= this->curr_size)) {
    cerr << "IndexedHashTable::member - i >= current size or < 0" << endl;
    exit(1);
  }

  ptr = IDX[i];
  tptr = ptr->item;

  return *tptr;
}

//
// IndexedHashTable<T>::member()
//
// Task:
//      return the ith element (sequentially).  Reference can be modified.
//      Use at your own peril!
//
// Conditions:
//      the index i must satisfy 0 <= i < curr_size
//

template <class T> T &IndexedHashTable<T>::member_ref(long i) {
  HEntry<T> *ptr;
  T *tptr;

  if ((i < 0) || (i >= this->curr_size)) {
    cerr << "IndexedHashTable::member - i >= current size or < 0" << endl;
    exit(1);
  }

  ptr = IDX[i];
  tptr = ptr->item;

  return *tptr;
}

//
// IndexedHashTable<T>::remove()
//
// Task:
//      removes G from the table, if it is there.
//

template <class T> void IndexedHashTable<T>::remove(const T &G) {
  long j;

  HEntry<T> *ptr;
  T *target, *tptr;

  // check whether G is in the table
  target = search(G);

  // if so, delete it
  if (target) {
    // compute j, its sequential index
    for (j = 0; j < this->curr_size; ++j) {
      ptr = IDX[j];
      tptr = ptr->item;
      if (tptr == target)
        break;
    }

    // remove element j
    remove_from(j);
  }
}

//
// IndexedHashTable<T>::remove_from()
//
// Task:
//      remove the item with sequential index i from the table
//
// Conditions:
//      the index i must satisfy 0 <= i < curr_size
//

template <class T> void IndexedHashTable<T>::remove_from(long i) {
  HEntry<T> *ptr, *pptr, *nptr, *temp;
  long j;

  if ((i < 0) || (i >= this->curr_size)) {
    cerr << "IndexedHashTable::remove_from - i >= current size or < 0" << endl;
    exit(1);
  }

  ptr = IDX[i];
  for (j = i; j < this->curr_size; ++j)
    IDX[j] = IDX[j + 1];
  IDX[j] = (HEntry<T> *)NULL;

  // find previous pointer in linked list
  j = hash_value(*(ptr->item)) % this->size;
  if (j < 0)
    j += this->size;
  temp = this->buckets[j];
  pptr = NULL;
  while (temp != ptr) {
    pptr = temp;
    temp = temp->next;
  }

  // delete
  nptr = ptr->next;
  if (pptr)
    pptr->next = nptr;
  else
    // G was the first element in the list - handle specially
    this->buckets[j] = nptr;

  delete ptr;

  --(this->curr_size);
}

//
// IndexedHashTable<T>::empty()
//
// Task:
//      delete all the elements in the table.  At the end, the number of
//      buckets is the same, but they will all be empty.
//

template <class T> void IndexedHashTable<T>::empty() {
  long i, end;

  end = this->curr_size;
  for (i = end - 1; i >= 0; --i)
    remove_from(i);
}

//
// IndexedHashTable<T>::hash()
//
// Task:
//      insert G into the hash table
//

template <class T> void IndexedHashTable<T>::hash(const T &G) {
  long i;

  HEntry<T> *ptr, *nptr, *newone;
  T *newT;

  // allocate new list element and a copy of G
  newone = new HEntry<T>;
  newT = new T;

  (*newT) = G;
  newone->item = newT;
  this->last_one = newT;

  // compute correct bucket and append G to the end of the linked list
  i = rem(hash_value(G), this->size);
  if (i < 0)
    i += this->size;
  ptr = this->buckets[i];
  if (ptr) {
    while (ptr) {
      nptr = ptr;
      ptr = nptr->next;
    }
    nptr->next = newone;
  } else
    this->buckets[i] = newone;

  // append the pointer to the list of pointers
  if (this->curr_size == allocated) {
    // allocated more storage for the sequential list

    HEntry<T> **tmp = new HEntry<T> *[allocated << 1];

    for (i = 0; i < allocated; ++i)
      tmp[i] = IDX[i];

    if (IDX)
      delete[] IDX;
    IDX = tmp;
    allocated <<= 1;
  }

  IDX[this->curr_size] = newone;
  ++(this->curr_size);
}

//
// IndexedHashTable<T>::output_style()
//
// Task:
//      set the style of output (0=list, 1=hash table)
//

template <class T> void IndexedHashTable<T>::output_style(int style) {
  if (style == 1)
    outstyle = 1;
  else
    outstyle = 0;
}

//
// IndexedHashTable<T>::read()
//
// Task:
//      read in a hash table from the istream in.
//
// Conditions:
//      input must be the number of buckets on one line, followed by the
//      number of elements to insert into the list on a new line, and finally
//      n instances of type T.
//

template <class T> void IndexedHashTable<T>::read(std::istream &in) {
  long i;
  long new_size, num;
  T new_item;

  in >> new_size;
  this->initialize((long)new_size);

  in >> num;
  for (i = 0; i < num; ++i) {
    in >> new_item;
    hash(new_item);
  }
}

//
// IndexedHashTable<T>::print()
//
// Task:
//      output a hash table to the ostream out.
//

template <class T> void IndexedHashTable<T>::print(std::ostream &out) const {
  long i;

  HEntry<T> *ptr;
  T *tptr;

  if (outstyle == 0) {
    // list style output
    for (i = 0; i < this->curr_size; ++i) {
      ptr = IDX[i];
      tptr = ptr->item;
      out << i << ": " << (*tptr) << "\n";
    }
  } else {
    // hash table style output
    for (i = 0; i < this->size; ++i) {
      if (this->buckets[i]) {
        out << "[" << i << "]";
        ptr = this->buckets[i];
        while (ptr) {
          tptr = ptr->item;
          out << " : " << (*tptr);
          ptr = ptr->next;
        }
        out << "\n";
      }
    }
  }
}

} // namespace ANTL
