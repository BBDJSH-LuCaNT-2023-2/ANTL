/**
 * @file qo_reduce.hpp
 * @author Michael Jacobson
 *
 * @brief Class defining abstract methods for ideal reduction
 */

#ifndef QO_REDUCE_H
#define QO_REDUCE_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/XGCD/xgcd.hpp>
#include <ANTL/Arithmetic/mul_exact.hpp>

#define RED_CANTOR 0
#define RED_FAST 1

NTL_CLIENT

template < class T> class QuadraticIdealBase;


//
// Class: qo_reduce<T>
//

template < class T > class qo_reduce
{
protected:
  T Delta;                         
  T hx;
  long genus;
  bool is_init;

public:
  qo_reduce() { is_init = false; };
  virtual ~qo_reduce() {};

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
  // Generic ideal reduction definition
  //
  virtual void reduce(QuadraticIdealBase<T> & A) = 0;
};

#endif // guard
