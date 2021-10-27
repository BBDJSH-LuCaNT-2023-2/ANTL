/**
 * @file qo_nudupl.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_square. 
 *  Computes a reduced ideal equivalent to the square of an ideal using the 
 *  NUDUPL algorithm, originally due to Shanks.
 */

#ifndef SQUARE_NUDUPL_H
#define SQUARE_NUDUPL_H

#include <ANTL/Quadratic/Square/SquareStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class SquareNudupl : public SquareStrategy<T> {
    using SquareStrategy<T>::Delta;
    using SquareStrategy<T>::hx;
    using SquareStrategy<T>::genus;
    using SquareStrategy<T>::is_init;

    protected:
      ZZ NC_BOUND;	  // termination bound for NUCOMP = floor(|D|^1/4)

    public:
      ~SquareNudupl() { };

      void init(const T & Din, const T & hin, long gin=0) {
        SquareStrategy<T>::init(Din,hin,gin);
      };

      // nudupl
      void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void SquareNudupl<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin);
template <> void SquareNudupl<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);

template <> void SquareNudupl<long>::init(const long & Din, const long & hin, long gin);
template <> void SquareNudupl<long>::square(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);

template <> void SquareNudupl<GF2EX>::square(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);

} //ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Square/SquareNudupl_impl.hpp"

#endif // guard
