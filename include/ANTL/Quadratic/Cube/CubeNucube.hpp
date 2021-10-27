/**
 * @file qo_nucube.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube. 
 *  Computes a reduced ideal equivalent to the cube of an ideal using the 
 *  NUCUBE algorithm.
 */

#ifndef CUBE_NUCUBE_H
#define CUBE_NUCUBE_H

#include <ANTL/Quadratic/Cube/CubeStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class CubeNucube : public CubeStrategy<T> {
    using CubeStrategy<T>::Delta;
    using CubeStrategy<T>::hx;
    using CubeStrategy<T>::genus;
    using CubeStrategy<T>::is_init;

    protected:
      ZZ SQRT_DELTA;   // = floor(SquareRoot(abs(Delta)))

    public:
      ~CubeNucube() { };

      void init(const T & Din, const T & hin, long gin=0) {
        CubeStrategy<T>::init(Din,hin,gin);
      };

      // nucube
      void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void CubeNucube<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin);
template <> void CubeNucube<ZZ>::cube(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);

template <> void CubeNucube<long>::init(const long & Din, const long & hin, long gin);
template <> void CubeNucube<long>::cube(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);

template <> void CubeNucube<GF2EX>::cube(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);

} //ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Cube/CubeNucube_impl.hpp"

#endif // guard
