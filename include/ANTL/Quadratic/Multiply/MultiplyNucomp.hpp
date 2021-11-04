/**
 * @file qo_nucomp.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_multiply. 
 *  Computes a reduced ideal equivalent to the product of two ideals using the 
 *  NUCOMP algorithm, originally due to Shanks.
 */

#ifndef MULTIPLY_NUCOMP_H
#define MULTIPLY_NUCOMP_H

#include <ANTL/Quadratic/Multiply/MultiplyStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT
using namespace ANTL;

namespace ANTL {

  template <class T> class MultiplyStrategy;
  template <class T> class QuadraticIdealBase;


  template < class T > class MultiplyNucomp : public MultiplyStrategy<T> {

    using MultiplyStrategy<T>::Delta;
    using MultiplyStrategy<T>::hx;
    using MultiplyStrategy<T>::genus;
    using MultiplyStrategy<T>::is_init;

    protected:
      ZZ NC_BOUND;    // termination bound for NUCOMP = floor(|D|^1/4)

    public:
      ~MultiplyNucomp() { };

    void init(const T & delta_in, const T & h_in, long g_in=0) {
      MultiplyStrategy<T>::init(delta_in,h_in,g_in);
    };

//     nucomp();
//     Task:
//          computes a reduced ideal equivalent to the product of two ideals
//          using the NUCOMP algorithm of Shanks.

     void multiply(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A, const QuadraticIdealBase<T> & B);
  };


//
// Declare specialized methods
//

template <> void MultiplyNucomp<ZZ>::init(const ZZ & delta_in, const ZZ & h_in, long g_in);
template <> void MultiplyNucomp<long>::init(const long & delta_in, const long & h_in, long g_in);

template <> void MultiplyNucomp<ZZ>::multiply(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A, const QuadraticIdealBase<ZZ> & B);

template <> void MultiplyNucomp<long>::multiply(QuadraticIdealBase<long> & C, const QuadraticIdealBase<long> & A, const QuadraticIdealBase<long> & B);

template <> void MultiplyNucomp<GF2EX>::multiply(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A, const QuadraticIdealBase<GF2EX> & B);

}//ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/Multiply/MultiplyNucomp_impl.hpp"

#endif // guard
