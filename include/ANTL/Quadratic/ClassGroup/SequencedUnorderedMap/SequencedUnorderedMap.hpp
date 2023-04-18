#ifndef ANTL_SEQUENCED_UNORDERED_MAP_H
#define ANTL_SEQUENCED_UNORDERED_MAP_H

NTL_CLIENT;

#include <unordered_map>
#include <list>
#include <queue>
#include <vector>

#include <ANTL/HashTable/HashEntryInt.hpp>

template <> struct std::hash<ANTL::HashEntryInt<long, ZZ>> {
  std::size_t operator()(ANTL::HashEntryInt<long, ZZ> const &hash_entry_int) const noexcept {
    std::size_t h1 = std::hash<long>{}(hash_entry_int.get_a());
    std::size_t h2 = std::hash<long>{}(hash_entry_int.get_b());
    return h1 ^ (h2 << 1); // or use boost::hash_combine
  }
};


template <class T> class SequencedUnorderedMap {
private:
  std::unordered_map<ANTL::HashEntryInt<T, ZZ>, std::queue<ANTL::HashEntryInt<T, ZZ>>> map;
  std::vector<ANTL::HashEntryInt<T, ZZ>> sequence;

public:
  SequencedUnorderedMap();
  ~SequencedUnorderedMap();

  void hash(ANTL::HashEntryInt<T, ZZ> hash_entry_int);
//   bool search(ANTL::HashEntryInt<T, ZZ> hash_entry_int);
  ANTL::HashEntryInt<T, ZZ> operator[](int i);
  ANTL::HashEntryInt<T, ZZ> * search(ANTL::HashEntryInt<T, ZZ> hash_entry_int);
  long no_of_elements();
  void remove_from(int i);
};

template <class T> SequencedUnorderedMap<T>::SequencedUnorderedMap() {}

template <class T> SequencedUnorderedMap<T>::~SequencedUnorderedMap() {}

template <class T> void SequencedUnorderedMap<T>::hash(ANTL::HashEntryInt<T, ZZ> hash_entry_int) {
  if(map.find(hash_entry_int) == map.end()) {
    map[hash_entry_int];
    map.at(hash_entry_int).push(hash_entry_int);
  }
  else {
    map.at(hash_entry_int).push(hash_entry_int);
  }
  sequence.push_back(hash_entry_int);
}

// template <class T> bool SequencedUnorderedMap<T>::search(ANTL::HashEntryInt<T, ZZ> hash_entry_int) {
//   if(map.find(hash_entry_int) != map.end()) {
//     return true;
//   }
//   return false;
// }

template <class T> ANTL::HashEntryInt<T, ZZ> SequencedUnorderedMap<T>::operator[](int i) {
  auto itr = sequence.begin();
  while(i > 0) {
    ++itr;
    --i;
  }
  return *itr;
}

template <class T> ANTL::HashEntryInt<T, ZZ> * SequencedUnorderedMap<T>::search(ANTL::HashEntryInt<T, ZZ> hash_entry_int) {
  if(map.find(hash_entry_int) == map.end()) {
    return (ANTL::HashEntryInt<T, ZZ> *)NULL;
  }
  else {
    return &map.at(hash_entry_int).front();
  }
}

template <class T> long SequencedUnorderedMap<T>::no_of_elements() {
  return sequence.size();
}

template <class T> void SequencedUnorderedMap<T>::remove_from(int i) {
  assert(i < this->no_of_elements());
    auto itr = sequence.begin();
  while(i > 0) {
    ++itr;
    --i;
  }

  map.at(*itr).pop();
  if(map.at(*itr).size() <= 0) {
    map.erase(*itr);
  }

  sequence.erase(itr);
}

#endif
