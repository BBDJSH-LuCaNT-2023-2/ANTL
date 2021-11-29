/**
 * @file qo_multiply_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_multiply.  Implements generic ideal multiplication.
 */

#ifndef MULTIPLY_PLAIN_H
#define MULTIPLY_PLAIN_H

#include <ANTL/Quadratic/Multiply/MultiplyStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>


NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class MultiplyStrategy;
  template <class T> class QuadraticIdealBase;

template < class T > class MultiplyPlain : public MultiplyStrategy<T>
{
  using MultiplyStrategy<T>::Delta;
  using MultiplyStrategy<T>::hx;
  using MultiplyStrategy<T>::genus;
  using MultiplyStrategy<T>::is_init;

  public:
  ~MultiplyPlain() { };

  //
  // generic ideal multiplication
  //
  void multiply(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A, const QuadraticIdealBase<T> & B);
};


//
// Declare specialized methods
//

template <> void MultiplyPlain<ZZ>::multiply(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A, const QuadraticIdealBase<ZZ> & B);

template <> void MultiplyPlain<long>::multiply(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A, const QuadraticIdealBase<long> & B);

template <> void MultiplyPlain<GF2EX>::multiply(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A, const QuadraticIdealBase<GF2EX> & B);


// Unspecialized template definitions.
#include "../src/Quadratic/Multiply/MultiplyPlain_impl.hpp"

} // ANTL
#endif // guard
