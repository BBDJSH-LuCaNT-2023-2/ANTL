/**
 * @file QuadraticOrder.hpp
 * @author Michael Jacobson
 * @brief order in a quadratic number field or hyperelliptic function field
 */

#ifndef ANTL_QUADRATIC_INF_ELEMENT_H
#define ANTL_QUADRATIC_INF_ELEMENT_H

#include <ANTL/Quadratic/QuadraticNumber.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

namespace ANTL {

using namespace ANTL;

  template < class T > class QuadraticNumber;
  template < class T > class QuadraticIdealBase;

  /**
   * @brief Element of a Quadratic Infrastructure
   * @remarks Uses a QuadraticIdealBase<T> and QuadraticNumber<T> to represent the
   *          and its distance along the Infrastructure respectively.
   */

template < class T > class QuadraticInfElement {
  private:
  QuadraticIdealBase<T> qib;
  RR Distance;

  T Delta;
  T FloorRootDelta;


  public:
  QuadraticInfElement();
  QuadraticInfElement(QuadraticOrder<T> & quad_o);

  ~QuadraticInfElement();

  QuadraticIdealBase<T> get_qib();
  RR get_distance();

  void baby_step();
  void giant_step(const QuadraticIdealBase<T> & quad_ib);
};

}
// Unspecialized template definitions.

template <> void QuadraticInfElement<ZZ>::baby_step();
template <> void QuadraticInfElement<long>::baby_step();

#include "../../../src/Quadratic/QuadraticInfElement_impl.hpp"

#endif // guard
