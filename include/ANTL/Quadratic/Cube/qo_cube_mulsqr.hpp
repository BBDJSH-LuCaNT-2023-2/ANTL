/**
 * @file qo_cube_mulsqr.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_cube.  
 *  Implements cubing by a square followed by multiplication
 */

#ifndef QO_CUBE_MULSQR_H
#define QO_CUBE_MULSQR_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Cube/qo_cube.hpp>

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
  void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A); //{
//    QuadraticIdealBase<T> CC(*A.get_QO());
//    sqr(CC,A);
//    mul(C,CC,A);
//  };
};

#include <../src/Quadratic/Cube/qo_cube_mulsqr_impl.hpp>

#endif // guard
