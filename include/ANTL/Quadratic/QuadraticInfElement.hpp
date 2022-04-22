/**
 * @file QuadraticOrder.hpp
 * @author Michael Jacobson
 * @brief order in a quadratic number field or hyperelliptic function field
 */

#ifndef ANTL_QUADRATIC_INF_ELEMENT_H
#define ANTL_QUADRATIC_INF_ELEMENT_H

#include <ANTL/HashTable/HashEntryReal.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>
#include <ANTL/Quadratic/QuadraticDistance.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/QuadraticNumber.hpp>

namespace ANTL {

using namespace ANTL;

template <class T> class QuadraticNumber;
template <class T> class QuadraticIdealBase;
template <class T, class S> class QuadraticInfElement;

// Forward declarations of friend functions

template <class T, class S>
void mul(QuadraticInfElement<T, S> qie_a, QuadraticInfElement<T, S> qie_b);
template <class T, class S>
void mul(QuadraticInfElement<T, S> qie_a, QuadraticInfElement<T, S> qie_b,
         QuadraticInfElement<T, S> qie_c);

template <class T, class S>
void sqr(QuadraticInfElement<T, S> qie_a, QuadraticInfElement<T, S> qie_b);
template <class T, class S>
void sqr(QuadraticInfElement<T, S> qie_a, QuadraticInfElement<T, S> qie_b,
         QuadraticInfElement<T, S> qie_c);

template <class T, class S>
void conjugate(QuadraticInfElement<T, S> &qie_a,
               QuadraticInfElement<T, S> const &qie_b);

template <class T, class S>
void nuclose(QuadraticInfElement<T, S> &C, const ZZ &n);

template <class S>
void update_distance_add(S &d_new, S const &distance_1, S const &distance_2);

template <class S>
void update_distance_subtract(S &d_new, S const &distance_1,
                              S const &distance_2);

template <class S>
void update_distance_multiply(S &d_new, S const &distance_1,
                              S const &distance_2);

template <class S> void update_distance_invert(S &d_new, S const &distance_1);

template <class S>
void update_distance_negate(S &distance_new, S const &distance_1);

/**
 * @brief Element of a Quadratic Infrastructure
 * @remarks Uses a QuadraticIdealBase<T> and RR to represent the
 *          and its distance along the Infrastructure respectively.
 */

template <class T, class S> class QuadraticInfElement {
private:
  QuadraticIdealBase<T> qib;
  S Distance;

  T Delta;
  T FloorRootDelta;

public:
  QuadraticInfElement();
  QuadraticInfElement(QuadraticOrder<T> &quad_o);

  ~QuadraticInfElement();

  QuadraticIdealBase<T> get_qib() const;
  S get_distance() const;

  // Friend functions for arithmetic
  friend void mul<T, S>(QuadraticInfElement<T, S> qie_a,
                        QuadraticInfElement<T, S> qie_b);
  friend void mul<T, S>(QuadraticInfElement<T, S> qie_a,
                        QuadraticInfElement<T, S> qie_b,
                        QuadraticInfElement<T, S> qie_c);

  friend void sqr<T, S>(QuadraticInfElement<T, S> qie_a,
                        QuadraticInfElement<T, S> qie_b);
  friend void sqr<T, S>(QuadraticInfElement<T, S> qie_a,
                        QuadraticInfElement<T, S> qie_b,
                        QuadraticInfElement<T, S> qie_c);

  // Friend function for Updating Distance
  friend void update_distance_add<S>(S &d_new, S const &distance_1,
                                     S const &distance_2);

  friend void update_distance_subtract<S>(S &d_new, S const &distance_1,
                                          S const &distance_2);

  friend void update_distance_multiply<S>(S &d_new, S const &distance_1,
                                          S const &distance_2);

  friend void update_distance_invert<S>(S &d_new, S const &distance_1);

  friend void update_distance_negate<S>(S &distance_new, S const &distance_1);

  friend void nuclose<T, S>(QuadraticInfElement<T, S> &C, const ZZ &n);

  // Infrastructure methods
  void baby_step();
  void giant_step(const QuadraticInfElement<T, S> &quad_ib);

  // START: TEMPORARY SECTION FOR REGULATORLENSTRADATA METHODS
  void adjust(const ZZ &a);
  void adjust(const S &a);
  void assign(const HashEntryReal<T, S> &her_a);
  void assign_one();

  friend void conjugate<T, S>(QuadraticInfElement<T, S> &qie_a,
                              QuadraticInfElement<T, S> const &qie_b);

  QuadraticInfElement<T, S> conjugate();

  ZZ eval();
  HashEntryReal<T, S> hash_real() const;
  bool is_one() const;
  void inverse_rho();

  S get_baby_steps(IndexedHashTable<HashEntryReal<T, S>> &prin_list,
                   const ZZ &B, const QuadraticInfElement<T, S> &A);

  S get_baby_steps(IndexedHashTable<HashEntryReal<T, S>> &prin_list,
                   const ZZ &B, const QuadraticInfElement<T, S> &A, long l,
                   long &M);

  //   qo_distance<T>
  //   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
  //                  const ZZ &B, const qi_pair<T> &A);
  //   qo_distance<T>
  //   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
  //                  const ZZ &B, const qi_pair<T> &A, long l, long &M);

  // FINISH: TEMPORARY SECTION FOR REGULATORLENSTRADATA METHODS
};
// Forward declarations of specialized template definitions.
//
// Special Distance Arithmetic Example
// template <>
// void update_distance_multiply<double>(double const &d_new,
//                                       double const &distance_1,
//                                       double const &distance_2);

} // namespace ANTL

#include "../../../src/Quadratic/QuadraticInfElement_impl.hpp"

#endif // guard
