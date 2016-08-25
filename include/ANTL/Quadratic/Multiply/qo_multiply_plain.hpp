/**
 * @file qo_multiply_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_multiply.  Implements generic ideal multiplication.
 */

#ifndef QO_MULTIPLY_PLAIN_H
#define QO_MULTIPLY_PLAIN_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Multiply/qo_multiply.hpp>

NTL_CLIENT


//
// Class: qo_multiply_plain<T>
//
template < class T > class qo_multiply_plain : public qo_multiply<T>
{
  using qo_multiply<T>::Delta;
  using qo_multiply<T>::hx;
  using qo_multiply<T>::genus;
  using qo_multiply<T>::is_init;

 public:
  ~qo_multiply_plain() { };

  //
  // generic ideal multiplication
  //
  void multiply(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A, const QuadraticIdealBase<T> & B);
};


//
// Declare specialized methods
//

template <> void qo_multiply_plain<ZZ>::multiply(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A, const QuadraticIdealBase<ZZ> & B);

template <> void qo_multiply_plain<GF2EX>::multiply(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A, const QuadraticIdealBase<GF2EX> & B);


// Unspecialized template definitions.
#include "../src/Quadratic/Multiply/qo_multiply_plain_impl.hpp"

#endif // guard
