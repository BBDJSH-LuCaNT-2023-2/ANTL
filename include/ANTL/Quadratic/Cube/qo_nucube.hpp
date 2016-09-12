/**
 * @file qo_nucube.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube. 
 *  Computes a reduced ideal equivalent to the cube of an ideal using the 
 *  NUCUBE algorithm.
 */

#ifndef QO_NUCUBE_H
#define QO_NUCUBE_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Cube/qo_cube.hpp>

NTL_CLIENT


//
// Class: qo_nucube<T>
//

template < class T > class qo_nucube : public qo_cube<T>
{
  using qo_cube<T>::Delta;
  using qo_cube<T>::hx;
  using qo_cube<T>::genus;
  using qo_cube<T>::is_init;

protected:
  ZZ SQRT_DELTA;   // = floor(SquareRoot(abs(Delta)))

 public:
  ~qo_nucube() { };

  void init(const T & Din, const T & hin, long gin=0) {
    qo_cube<T>::init(Din,hin,gin);
  };

  //
  // nucube
  //
  void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_nucube<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin);
template <> void qo_nucube<ZZ>::cube(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);

template <> void qo_nucube<long>::init(const long & Din, const long & hin, long gin);
template <> void qo_nucube<long>::cube(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A);

template <> void qo_nucube<GF2EX>::cube(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);


// Unspecialized template definitions.
#include "../src/Quadratic/Cube/qo_nucube_impl.hpp"

#endif // guard
