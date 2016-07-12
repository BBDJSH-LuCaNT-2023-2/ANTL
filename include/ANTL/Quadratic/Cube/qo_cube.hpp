/**
 * @file qo_cube.hpp
 * @author Michael Jacobson
 * @brief  Class defining abstract methods for ideal cubing.
 */

#ifndef QO_CUBE_H
#define QO_CUBE_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>

#include <Quadratic/QuadraticIdealBase.hpp>
#include <XGCD/xgcd.hpp>
#include <Arithmetic/mul_exact.hpp>

#define CUBE_CANTOR 0
#define CUBE_NUCOMP 1
#define CUBE_MULSQR 2
#define CUBE_EXPLICIT 3

NTL_CLIENT

template < class T> class QuadraticIdealBase;

//
// Class: qo_cube<T>
//
template < class T > class qo_cube
{
protected:
  T Delta;                         
  T hx;
  long genus;
  bool is_init;

public:
  qo_cube() { is_init=false; };
  virtual ~qo_cube() {};
 
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
  // Generic ideal cubing definition
  //
  virtual void cube(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A) = 0;
};

#endif // guard
