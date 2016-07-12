/**
 * @file qo_nucomp.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_multiply. 
 *  Computes a reduced ideal equivalent to the product of two ideals using the 
 *  NUCOMP algorithm, originally due to Shanks.
 */

#ifndef QO_NUCOMP_H
#define QO_NUCOMP_H

#include <Quadratic/QuadraticIdealBase.hpp>
#include <Quadratic/Multiply/qo_multiply.hpp>

NTL_CLIENT

//
// Class: qo_nucomp<T>
//

template < class T > class qo_nucomp : public qo_multiply<T>
{
  using qo_multiply<T>::Delta;
  using qo_multiply<T>::hx;
  using qo_multiply<T>::genus;
  using qo_multiply<T>::is_init;

protected:
  ZZ NC_BOUND;	  // termination bound for NUCOMP = floor(|D|^1/4)

 public:
  ~qo_nucomp() { };

  void init(T Din, T hin, long gin=0) {
    qo_multiply<T>::init(Din,hin,gin);
  };

  //
  // nucomp();
  //
  // Task:
  //      computes a reduced ideal equivalent to the product of two ideals
  //      using the NUCOMP algorithm of Shanks.
  //
  void multiply(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A, const QuadraticIdealBase<T> & B);
};


//
// Declare specialized methods
//

template <> void qo_nucomp<ZZ>::init(ZZ Din, ZZ hin, long gin);
template <> void qo_nucomp<ZZ>::multiply(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A, const QuadraticIdealBase<ZZ> & B);
template <> void qo_nucomp<GF2EX>::multiply(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A, const QuadraticIdealBase<GF2EX> & B);


// Unspecialized template definitions.
#include "../src/Quadratic/Multiply/qo_nucomp_impl.hpp"

#endif // guard
