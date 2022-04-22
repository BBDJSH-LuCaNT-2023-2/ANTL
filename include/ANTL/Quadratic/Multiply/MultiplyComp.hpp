/**
 * @file qo_multiply_plain.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_multiply.  Implements generic ideal multiplication.
 */

#ifndef MULTIPLY_COMP_H
#define MULTIPLY_COMP_H

#include <ANTL/Quadratic/Multiply/MultiplyStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>


NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class MultiplyStrategy;
  template <class T> class QuadraticIdealBase;
  template <class T> class QuadraticNumber;

template < class T > class MultiplyComp : public MultiplyStrategy<T>
{
  using MultiplyStrategy<T>::Delta;
  using MultiplyStrategy<T>::hx;
  using MultiplyStrategy<T>::genus;
  using MultiplyStrategy<T>::is_init;

  protected:
      QuadraticNumber<T> * RelativeGenerator;

  public:
  ~MultiplyComp() { };

    QuadraticNumber<T> * get_RelativeGenerator() {
      return RelativeGenerator;
    }

    void set_RelativeGenerator(QuadraticNumber<T> & QN) {
        RelativeGenerator = &QN;
    }

  //
  // generic ideal multiplication
  //
  void multiply(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A, const QuadraticIdealBase<T> & B);
};


//
// Declare specialized methods
//

template <> void MultiplyComp<ZZ>::multiply(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A, const QuadraticIdealBase<ZZ> & B);

template <> void MultiplyComp<long>::multiply(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A, const QuadraticIdealBase<long> & B);

template <> void MultiplyComp<GF2EX>::multiply(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A, const QuadraticIdealBase<GF2EX> & B);


// Unspecialized template definitions.
#include "../src/Quadratic/Multiply/MultiplyComp_impl.hpp"

} // ANTL
#endif // guard
