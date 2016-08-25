/**
 * @file qo_nudupl.hpp
 * @author Michael Jacobson
 * @brief Concrete class extending qo_square. 
 *  Computes a reduced ideal equivalent to the square of an ideal using the 
 *  NUDUPL algorithm, originally due to Shanks.
 */

#ifndef QO_NUDUPL_H
#define QO_NUDUPL_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Square/qo_square.hpp>

NTL_CLIENT


//
// Class: qo_nudupl<T>
//

template < class T > class qo_nudupl : public qo_square<T>
{
  using qo_square<T>::Delta;
  using qo_square<T>::hx;
  using qo_square<T>::genus;
  using qo_square<T>::is_init;

protected:
  ZZ NC_BOUND;	  // termination bound for NUCOMP = floor(|D|^1/4)

 public:
  ~qo_nudupl() { };

  void init(T Din, T hin, long gin=0) {
    qo_square<T>::init(Din,hin,gin);
  };

  //
  // nudupl
  //
  void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A);
};


//
// Declare specialized methods
//

template <> void qo_nudupl<ZZ>::init(ZZ Din, ZZ hin, long gin);
template <> void qo_nudupl<ZZ>::square(QuadraticIdealBase<ZZ> & C, const QuadraticIdealBase<ZZ> & A);
template <> void qo_nudupl<GF2EX>::square(QuadraticIdealBase<GF2EX> & C, const QuadraticIdealBase<GF2EX> & A);



// Unspecialized template definitions.
#include "../src/Quadratic/Square/qo_nudupl_impl.hpp"

#endif // guard
