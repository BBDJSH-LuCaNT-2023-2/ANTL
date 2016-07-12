/**
 * @file qo_cube_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube.  Implements generic ideal cubing.
 */

#ifndef QO_CUBE_PLAIN_H
#define QO_CUBE_PLAIN_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Cube/qo_cube.hpp>

NTL_CLIENT


//
// Class: qo_cube_plain<T>
//
template < class T > class qo_cube_plain : public qo_cube<T>
{
  using qo_cube<T>::Delta;
  using qo_cube<T>::hx;
  using qo_cube<T>::genus;
  using qo_cube<T>::is_init;

 public:
  ~qo_cube_plain() { };

  //
  // generic ideal cubing
  //
  void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_cube_plain<ZZ>::cube(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);
template <> void qo_cube_plain<GF2EX>::cube(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);


// Unspecialized template definitions.
#include "../src/Quadratic/Cube/qo_cube_plain_impl.hpp"

#endif // guard
