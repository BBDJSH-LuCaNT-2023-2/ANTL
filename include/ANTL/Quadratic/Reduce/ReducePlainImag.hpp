/**
 * @file qo_reduce_plain_imag.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implement's Cantor's reduction algorithm.
 */

#ifndef REDUCE_PLAIN_IMAG_H
#define REDUCE_PLAIN_IMAG_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReducePlainImag : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReducePlainImag() {};

      void reduce(QuadraticIdealBase<T> & A);
  };


//
// Declare specialized methods
//

template <> void ReducePlainImag<ZZ>::reduce(ANTL::QuadraticIdealBase<ZZ> & A);
template <> void ReducePlainImag<long>::reduce(ANTL::QuadraticIdealBase<long> & A);
template <> void ReducePlainImag<GF2EX>::reduce(ANTL::QuadraticIdealBase<GF2EX> & A);

} // ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReducePlainImag_impl.hpp"

#endif // guard
