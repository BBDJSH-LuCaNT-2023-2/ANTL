/**
 * @file qo_cube_mulsqr.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube.  
 *  Implements cubing by a square followed by multiplication
 */

#ifndef QO_CUBE_MULSQR_H
#define QO_CUBE_MULSQR_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Cube/qo_cube.hpp>

NTL_CLIENT


//
// Class: qo_cube_mulsqr<T>
//

template < class T > class qo_cube_mulsqr : public qo_cube<T>
{
  using qo_cube<T>::Delta;
  using qo_cube<T>::hx;
  using qo_cube<T>::genus;
  using qo_cube<T>::is_init;

 public:
  ~qo_cube_mulsqr() { };

  //
  // cubing via square and multiply
  //
  void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A) {
    QuadraticIdealBase<T> CC;
    square(CC,A);
    multiply(C,CC,A);
  };
};

#endif // guard
