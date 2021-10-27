/**
 * @file qo_reduce_fast_allrecs.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements and adaptation of Sawilla's
 * algorithm with all extra recurrences computed during the partial EEA.
 */

#ifndef REDUCE_FAST_ALLRECS_H
#define REDUCE_FAST_ALLRECS_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReduceFastAllRecs : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReduceFastAllRecs() {};

      void init(T Din, T hin, long gin=0) {
        ReduceStrategy<T>::init(Din,hin,gin);
      };

      void reduce(QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void ReduceFastAllRecs<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);

}//ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReduceFastAllRecs_impl.hpp"

#endif // guard
