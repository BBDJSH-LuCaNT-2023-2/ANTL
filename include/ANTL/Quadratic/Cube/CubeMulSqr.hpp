/**
 * @file qo_cube_mulsqr.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube.  
 *  Implements cubing by a square followed by multiplication
 */

#ifndef CUBE_MULSQR_H
#define CUBE_MULSQR_H

#include <ANTL/Quadratic/Cube/CubeStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class CubeMulSqr : public CubeStrategy<T> {
    using CubeStrategy<T>::Delta;
    using CubeStrategy<T>::hx;
    using CubeStrategy<T>::genus;
    using CubeStrategy<T>::is_init;

    public:
      ~CubeMulSqr() { };

      // Cubing via square and multiply
      void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
        //{
        // QuadraticIdealBase<T> CC(*A.get_QO());
        // sqr(CC,A);
        // mul(C,CC,A);
        // };
  };

} //ANTL
#include <../src/Quadratic/Cube/CubeMulSqr_impl.hpp>

#endif // guard
