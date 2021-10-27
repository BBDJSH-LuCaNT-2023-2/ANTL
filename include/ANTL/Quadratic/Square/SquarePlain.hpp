/**
 * @file qo_square_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_square.  Implements generic ideal squaring.
 */

#ifndef SQUARE_PLAIN_H
#define SQUARE_PLAIN_H

#include <ANTL/Quadratic/Square/SquareStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class SquarePlain : public SquareStrategy<T> {
    using SquareStrategy<T>::Delta;
    using SquareStrategy<T>::hx;
    using SquareStrategy<T>::genus;
    using SquareStrategy<T>::is_init;

    public:
      ~SquarePlain() { };

      // generic ideal squaring
      void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void SquarePlain<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);
template <> void SquarePlain<long>::square(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);
template <> void SquarePlain<GF2EX>::square(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);

} //ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Square/SquarePlain_impl.hpp"

#endif // guard
