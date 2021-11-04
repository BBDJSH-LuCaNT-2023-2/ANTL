/**
 * @file qo_reduce_fast_allrecs_fancy.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements and adaptation of Sawilla's
 * algorithm with all extra recurrences computed during the partial EEA.  
 * Coefficient c is computed directly from the recurrenes.
 */

#ifndef REDUCE_FAST_ALLRECS_FANCY_H
#define REDUCE_FAST_ALLRECS_FANCY_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReduceFastAllRecsFancy : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReduceFastAllRecsFancy() {};

      void init(T delta_in, T h_in, long g_in=0) {
        ReduceStrategy<T>::init(delta_in,h_in,g_in);
      };

      void reduce(QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void ReduceFastAllRecsFancy<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);

}//ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReduceFastAllRecsFancy_impl.hpp"

#endif // guard
