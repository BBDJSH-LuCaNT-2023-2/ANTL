#ifndef REDUCE_PLAIN_REAL_OPT_H
#define REDUCE_PLAIN_REAL_OPT_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReducePlainRealOpt : public ReduceStrategy<T> {
    using ReduceStrategy<T>::RelativeGenerator;
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::FloorRootDelta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    public:
      ~ReducePlainRealOpt() {};

      void reduce(QuadraticIdealBase<T> & A);

      using ReduceStrategy<T>::get_RelativeGenerator;
      using ReduceStrategy<T>::set_RelativeGenerator;
  };


// Declaration of specialized methods
template <> void ReducePlainRealOpt<ZZ>::reduce(ANTL::QuadraticIdealBase<ZZ> & A);
template <> void ReducePlainRealOpt<long>::reduce(ANTL::QuadraticIdealBase<long> & A);
template <> void ReducePlainRealOpt<GF2EX>::reduce(ANTL::QuadraticIdealBase<GF2EX> & A);

} // ANTL

// Unspecialized template definitions.

// TODO: impl file
// #include "../src/Quadratic/Reduce/ReducePlainRealOpt_impl.hpp"

#endif // guard

