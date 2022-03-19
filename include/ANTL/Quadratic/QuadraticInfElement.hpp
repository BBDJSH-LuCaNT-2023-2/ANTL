/**
 * @file QuadraticOrder.hpp
 * @author Michael Jacobson
 * @brief order in a quadratic number field or hyperelliptic function field
 */

#ifndef ANTL_QUADRATIC_INF_ELEMENT_H
#define ANTL_QUADRATIC_INF_ELEMENT_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/QuadraticNumber.hpp>

namespace ANTL {

using namespace ANTL;

template <class T> class QuadraticNumber;
template <class T> class QuadraticIdealBase;

/**
 * @brief Element of a Quadratic Infrastructure
 * @remarks Uses a QuadraticIdealBase<T> and RR to represent the
 *          and its distance along the Infrastructure respectively.
 */

template <class T> class QuadraticInfElement {
private:
  QuadraticIdealBase<T> qib;
  RR Distance;

  T Delta;
  T FloorRootDelta;

public:
  QuadraticInfElement();
  QuadraticInfElement(QuadraticOrder<T> &quad_o);

  ~QuadraticInfElement();

  QuadraticIdealBase<T> get_qib();
  RR get_distance();

  void baby_step();
  void giant_step(QuadraticInfElement<T> &quad_ib);
  void giant_step(const QuadraticIdealBase<T> &quad_ib);

  // TEMPORARY SECTION FOR REGULATORLENSTRADATA METHODS
  void adjust(const ZZ &a);
  void assign();
  void assign_one();
  bool is_one();

  RR get_baby_steps(IndexedHashTable<HashEntryReal<T>> &prin_list, const ZZ &B,
                    const QuadraticInfElement<T> &A);

  RR get_baby_steps(IndexedHashTable<HashEntryReal<T>> &prin_list, const ZZ &B,
                    const QuadraticInfElement<T> &A, long l, long &M);

  //   qo_distance<T>
  //   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
  //                  const ZZ &B, const qi_pair<T> &A);
  //   qo_distance<T>
  //   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
  //                  const ZZ &B, const qi_pair<T> &A, long l, long &M);

  // TEMPORARY SECTION FOR REGULATORLENSTRADATA METHODS
};

} // namespace ANTL
// Unspecialized template definitions.

template <> void QuadraticInfElement<ZZ>::baby_step();
template <> void QuadraticInfElement<long>::baby_step();

#include "../../../src/Quadratic/QuadraticInfElement_impl.hpp"

#endif // guard
