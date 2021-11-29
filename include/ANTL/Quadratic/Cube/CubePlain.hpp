/**
 * @file qo_cube_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube.  Implements generic ideal cubing.
 */

#ifndef CUBE_PLAIN_H
#define CUBE_PLAIN_H

#include <ANTL/Quadratic/Cube/CubeStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class CubePlain : public CubeStrategy<T> {
    using CubeStrategy<T>::Delta;
    using CubeStrategy<T>::hx;
    using CubeStrategy<T>::genus;
    using CubeStrategy<T>::is_init;

    public:
      ~CubePlain() {};

      void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
    };

  // Declare specialized methods
  template <> void CubePlain<ZZ>::cube(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);
  template <> void CubePlain<long>::cube(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);
  template <> void CubePlain<GF2EX>::cube(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);

} //ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Cube/CubePlain_impl.hpp"

#endif // guard
