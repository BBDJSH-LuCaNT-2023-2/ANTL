/**
 * @file qo_square.hpp
 * @author Michael Jacobson
 * @brief Class defining abstract methods for ideal squaring.
 */

#ifndef QO_SQUARE_H
#define QO_SQUARE_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>

#include <Quadratic/QuadraticIdealBase.hpp>
#include <XGCD/xgcd.hpp>
#include <Arithmetic/mul_exact.hpp>

#define SQR_CANTOR 0
#define SQR_NUCOMP 1
#define SQR_EXPLICIT 2

NTL_CLIENT

template < class T> class QuadraticIdealBase;

//
// Class: qo_square<T>
//

template < class T > class qo_square
{
protected:
  T Delta;                         
  T hx;
  long genus;
  bool is_init;

public:
  qo_square() { is_init = false; };
  virtual ~qo_square() {};

  //
  // Initialize field invariants
  //
  void init(T Din, T hin, long gin=0) {
    if (is_init) {
      Delta.kill();
      hx.kill();
    }

    Delta = Din;
    hx = hin;
    genus = gin;
    is_init = true;
  };

  //
  // Generic ideal squaring definition
  //
  virtual void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A) = 0;
};

#endif // guard
