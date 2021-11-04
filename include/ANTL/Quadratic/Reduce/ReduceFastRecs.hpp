/**
 * @file qo_reduce_fast_recs.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements an adaptation of Sawilla's 
 * reduction algorithm with some recurrences computed during the partial EEA.
 */

#ifndef REDUCE_FAST_RECS_H
#define REDUCE_FAST_RECS_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>


NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReduceFastRecs : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReduceFastRecs() {};

      void init(T delta_in, T h_in, long g_in=0) {
        ReduceStrategy<T>::init(delta_in,h_in,g_in);
      };

      void reduce(QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void ReduceFastRecs<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);

}//ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReduceFastRecs_impl.hpp"

#endif // guard
