/**
 * @file qo_reduce_fast_alt_a.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements an adaptation of Sawilla's 
 * algorithm with some recurrences computed during the partial EEA.  
 * Uses alternative formula for a.
 */

#ifndef REDUCE_FAST_ALT_A_H
#define REDUCE_FAST_ALT_A_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReduceFastAltA : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReduceFastAltA() {};

      void init(T delta_in, T h_in, long g_in=0) {
        ReduceStrategy<T>::init(delta_in,h_in,g_in);
      };

      void reduce(QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void ReduceFastAltA<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);

}//ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReduceFastAltA_impl.hpp"

#endif // guard
