/**
 * @file qo_reduce_plain_real.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implement's Cantor's reduction algorithm.
 */

#ifndef REDUCE_PLAIN_REAL_H
#define REDUCE_PLAIN_REAL_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReducePlainReal : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReducePlainReal() {};

      void normalize(QuadraticIdealBase<T> & A);
      void reduce(QuadraticIdealBase<T> & A);
  };


// Declaration of specialized methods
template <> void ReducePlainReal<ZZ>::reduce(ANTL::QuadraticIdealBase<ZZ> & A);
template <> void ReducePlainReal<long>::reduce(ANTL::QuadraticIdealBase<long> & A);
template <> void ReducePlainReal<GF2EX>::reduce(ANTL::QuadraticIdealBase<GF2EX> & A);

template <> void ReducePlainReal<ZZ>::normalize(ANTL::QuadraticIdealBase<ZZ> & A);
template <> void ReducePlainReal<long>::normalize(ANTL::QuadraticIdealBase<long> & A);

} // ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReducePlainReal_impl.hpp"

#endif // guard
