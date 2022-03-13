#ifndef REGULATOR_LENSTRA_DATA_H
#define REGULATOR_LENSTRA_DATA_H

#include <ANTL/L_function/L_function.hpp>
#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>
#include <ANTL/common.hpp>

NTL_CLIENT
using namespace ANTL;

namespace ANTL {
template <class T> class RegulatorLenstraData {
private:
  QuadraticOrder<T> quadratic_order;
  T delta;

  L_function<T> l_function;

  bool parallel;

  const long OQvals_cnum[20] = {2269,   5741,   10427,  16183,  22901,
                                30631,  39209,  48731,  59063,  70237,
                                82223,  95009,  108571, 122921, 137983,
                                153817, 170341, 187631, 205589, 224261};

  ZZ max_memory{4000000000};

  indexed_hash_table < qo_hash_entry_real < T > >prin_list;

  long prinlist_l;	// distance between consecutive table entries

  long prinlist_M;	// max # of baby-steps between table entries

  ZZ prinlist_s;		// giant-step distance

public:
  QuadraticOrder<T> get_quadratic_order();

  T get_delta();

  L_function<T> get_l_function();

  bool get_parallel();

  long get_OQvals_cnum_entry(long index);

  ZZ get_max_memory();

};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class T>
QuadraticOrder<T> RegulatorLenstraData<T>::get_quadratic_order() {
  return quadratic_order;
}

template <class T> T RegulatorLenstraData<T>::get_delta() { return delta; }

template <class T> L_function<T> RegulatorLenstraData<T>::get_l_function() {
  return l_function;
}

template <class T>
long RegulatorLenstraData<T>::get_OQvals_cnum_entry(long index) {
  return OQvals_cnum[index];
}

template <class T> bool RegulatorLenstraData<T>::get_parallel() {
  return parallel;
}

template <class T> ZZ RegulatorLenstraData<T>::get_max_memory() {
  return max_memory;
}
} // namespace ANTL
#endif
