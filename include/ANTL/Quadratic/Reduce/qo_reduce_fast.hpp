/**
 * @file qo_reduce_fast.hpp
 * @author Michael Jacobson
 *
 * @brief Concrete class extending qo_reduce.  Implements Sawilla's reduction algorithm
 * and its adaptation to function fields.
 */

#ifndef QO_REDUCE_FAST_H
#define QO_REDUCE_FAST_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Reduce/qo_reduce.hpp>

NTL_CLIENT

//
// Class: qo_reduce_fast<T>
//

template < class T > class qo_reduce_fast : public qo_reduce<T>
{
  using qo_reduce<T>::Delta;
  using qo_reduce<T>::hx;
  using qo_reduce<T>::genus;
  using qo_reduce<T>::is_init;

protected:
  ZZ SQRT_DELTA;   // = floor(SquareRoot(abs(Delta)))

 public:
  ~qo_reduce_fast() {};

  void init(T Din, T hin, long gin=0) {
    qo_reduce<T>::init(Din,hin,gin);
  };

  void reduce(QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_reduce_fast<ZZ>::init(ZZ Din, ZZ hin, long gin);
template <> void qo_reduce_fast<ZZ>::reduce(QuadraticIdealBase<ZZ> & A);
template <> void qo_reduce_fast<GF2EX>::reduce(QuadraticIdealBase<GF2EX> & A);


// Unspecialized template definitions.
#include "../src/Quadratic/Reduce/qo_reduce_fast_impl.hpp"

#endif // guard
