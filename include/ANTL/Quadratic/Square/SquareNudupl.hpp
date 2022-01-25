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

  template <class T> class SquareStrategy;
  template <class T> class QuadraticIdealBase;
  template <class T> class QuadraticNumber;

  template <class T> class SquareNudupl : public SquareStrategy<T> {
    using SquareStrategy<T>::Delta;
    using SquareStrategy<T>::hx;
    using SquareStrategy<T>::genus;
    using SquareStrategy<T>::is_init;

    protected:
      QuadraticNumber<T> * RelativeGenerator;
      ZZ NC_BOUND;	  // termination bound for NUCOMP = floor(|D|^1/4)

    public:
      ~SquareNudupl() { };

      void init(const T & delta_in, const T & h_in, long g_in=0) {
        SquareStrategy<T>::init(delta_in,h_in,g_in);
      };

      QuadraticNumber<T> * get_RelativeGenerator() {
        return RelativeGenerator;
      }

      void set_RelativeGenerator(QuadraticNumber<T> & QN) {
        RelativeGenerator = &QN;
      }

      // nudupl
      void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void SquareNudupl<ZZ>::init(const ZZ & delta_in, const ZZ & h_in, long g_in);
template <> void SquareNudupl<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);

template <> void SquareNudupl<long>::init(const long & delta_in, const long & h_in, long g_in);
template <> void SquareNudupl<long>::square(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);

template <> void SquareNudupl<GF2EX>::square(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);

} //ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Square/SquareNudupl_impl.hpp"

#endif // guard
