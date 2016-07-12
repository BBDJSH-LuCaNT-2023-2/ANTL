/**
 * @file qo_square_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_square.  Implements generic ideal squaring.
 */

#ifndef QO_SQUARE_PLAIN_H
#define QO_SQUARE_PLAIN_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Square/qo_square.hpp>

NTL_CLIENT


//
// Class: qo_square_plain<T>
//
template < class T > class qo_square_plain : public qo_square<T>
{
  using qo_square<T>::Delta;
  using qo_square<T>::hx;
  using qo_square<T>::genus;
  using qo_square<T>::is_init;

 public:
  ~qo_square_plain() { };

  //
  // generic ideal squaring
  //
  void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_square_plain<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);
template <> void qo_square_plain<GF2EX>::square(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);


// Unspecialized template definitions.
#include "../src/Quadratic/Square/qo_square_plain_impl.hpp"

#endif // guard
