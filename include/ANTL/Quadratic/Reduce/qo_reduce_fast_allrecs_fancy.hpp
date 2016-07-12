/**
 * @file qo_reduce_fast_allrecs_fancy.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements and adaptation of Sawilla's
 * algorithm with all extra recurrences computed during the partial EEA.  
 * Coefficient c is computed directly from the recurrenes.
 */

#ifndef QO_REDUCE_FAST_ALLRECS_FANCY_H
#define QO_REDUCE_FAST_ALLRECS_FANCY_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Reduce/qo_reduce.hpp>

NTL_CLIENT


//
// Class: qo_reduce_fast_allrecs_fancy<T>
//

template < class T > class qo_reduce_fast_allrecs_fancy : public qo_reduce<T>
{
  using qo_reduce<T>::Delta;
  using qo_reduce<T>::hx;
  using qo_reduce<T>::genus;
  using qo_reduce<T>::is_init;

 public:
  ~qo_reduce_fast_allrecs_fancy() {};

  void init(T Din, T hin, long gin=0) {
    qo_reduce<T>::init(Din,hin,gin);
  };

  void reduce(QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_reduce_fast_allrecs_fancy<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);


// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/qo_reduce_fast_allrecs_fancy_impl.hpp"

#endif // guard
