/**
 * @file qo_reduce_fast_alt_a.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements an adaptation of Sawilla's 
 * algorithm with some recurrences computed during the partial EEA.  
 * Uses alternative formula for a.
 */

#ifndef QO_REDUCE_FAST_ALT_A_H
#define QO_REDUCE_FAST_ALT_A_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Reduce/qo_reduce.hpp>

NTL_CLIENT


//
// Class: qo_reduce_fast_alt_a<T>
//

template < class T > class qo_reduce_fast_alt_a : public qo_reduce<T>
{
  using qo_reduce<T>::Delta;
  using qo_reduce<T>::hx;
  using qo_reduce<T>::genus;
  using qo_reduce<T>::is_init;

 public:
  ~qo_reduce_fast_alt_a() {};

  void init(T Din, T hin, long gin=0) {
    qo_reduce<T>::init(Din,hin,gin);
  };

  void reduce(QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_reduce_fast_alt_a<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);


// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/qo_reduce_fast_alt_a_impl.hpp"

#endif // guard
