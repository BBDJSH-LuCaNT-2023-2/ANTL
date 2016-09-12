/**
 * @file qo_reduce_plain_imag.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implement's Cantor's reduction algorithm.
 */

#ifndef QO_REDUCE_PLAIN_IMAG_H
#define QO_REDUCE_PLAIN_IMAG_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Reduce/qo_reduce.hpp>

NTL_CLIENT


//
// Class: qo_reduce_plain_imag<T>
//

template < class T > class qo_reduce_plain_imag : public qo_reduce<T>
{
  using qo_reduce<T>::Delta;
  using qo_reduce<T>::hx;
  using qo_reduce<T>::genus;
  using qo_reduce<T>::is_init;

 public:
  ~qo_reduce_plain_imag() {};

  void reduce(QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_reduce_plain_imag<ZZ>::reduce(QuadraticIdealBase<ZZ> & A);
template <> void qo_reduce_plain_imag<long>::reduce(QuadraticIdealBase<long> & A);
template <> void qo_reduce_plain_imag<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);



// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/qo_reduce_plain_imag_impl.hpp"

#endif // guard
