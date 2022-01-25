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

  template <class T> class CubeStrategy;
  template <class T> class QuadraticIdealBase;
  template <class T> class QuadraticNumber;

  template <class T> class CubeNucube : public CubeStrategy<T> {
    using CubeStrategy<T>::Delta;
    using CubeStrategy<T>::hx;
    using CubeStrategy<T>::genus;
    using CubeStrategy<T>::is_init;

    protected:
      QuadraticNumber<T> * RelativeGenerator;
      ZZ sqrt_delta;   // = floor(SquareRoot(abs(Delta)))

    public:
      ~CubeNucube() { };

      void init(const T & delta_in, const T & h_in, long g_in=0) {
        CubeStrategy<T>::init(delta_in,h_in,g_in);
      };

      QuadraticNumber<T> * get_RelativeGenerator() {
        return RelativeGenerator;
      }

      void set_RelativeGenerator(QuadraticNumber<T> & QN) {
        RelativeGenerator = &QN;
      }

      // nucube
      void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
  };

// Declare specialized methods
template <> void CubeNucube<ZZ>::init(const ZZ & delta_in, const ZZ & h_in, long g_in);
template <> void CubeNucube<ZZ>::cube(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);

template <> void CubeNucube<long>::init(const long & delta_in, const long & h_in, long g_in);
template <> void CubeNucube<long>::cube(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);

template <> void CubeNucube<GF2EX>::cube(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);

} //ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Cube/CubeNucube_impl.hpp"

#endif // guard
