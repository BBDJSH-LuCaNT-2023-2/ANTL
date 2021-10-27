/**
 * @file qo_reduce_fast.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements Sawilla's reduction algorithm
 * and its adaptation to function fields.
 */

#ifndef REDUCE_FAST_H
#define REDUCE_FAST_H

#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReduceFast : public ReduceStrategy<T> {
    using ReduceStrategy<T>::Delta;
    using ReduceStrategy<T>::hx;
    using ReduceStrategy<T>::genus;
    using ReduceStrategy<T>::is_init;

    protected:
      ZZ SQRT_DELTA;   // = floor(SquareRoot(abs(Delta)))

    public:
      ~ReduceFast() {};

      void init(const T & Din, const T & hin, long gin=0) {
        ReduceStrategy<T>::init(Din,hin,gin);
      }

      void reduce(QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void ReduceFast<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin);
template <> void ReduceFast<ZZ>::reduce(QuadraticIdealBase<ZZ> & A);

template <> void ReduceFast<long>::init(const long & Din, const long & hin, long gin);
template <> void ReduceFast<long>::reduce(QuadraticIdealBase<long> & A);

template <> void ReduceFast<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);

} //ANTL
// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/ReduceFast_impl.hpp"

#endif // guard
