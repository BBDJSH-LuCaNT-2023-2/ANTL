/**
 * @file qo_reduce_plain.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implement's Cantor's reduction algorithm.
 */

#ifndef QO_REDUCE_PLAIN_H
#define QO_REDUCE_PLAIN_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Reduce/qo_reduce.hpp>

NTL_CLIENT


//
// Class: qo_reduce_plain<T>
//

template < class T > class qo_reduce_plain : public qo_reduce<T>
{
  using qo_reduce<T>::Delta;
  using qo_reduce<T>::hx;
  using qo_reduce<T>::genus;
  using qo_reduce<T>::is_init;

 public:
  ~qo_reduce_plain() {};

  void reduce(QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_reduce_plain<ZZ>::reduce(QuadraticIdealBase<ZZ> & A);
template <> void qo_reduce_plain<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);



// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/qo_reduce_plain_impl.hpp"

#endif // guard
