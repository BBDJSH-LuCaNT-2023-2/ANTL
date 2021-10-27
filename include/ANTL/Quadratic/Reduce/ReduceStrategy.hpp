/**
 * @file qo_reduce.hpp
 * @author Michael Jacobson
 *
 * @brief Class defining abstract methods for ideal reduction
 */

#ifndef REDUCE_STRATEGY_H
#define REDUCE_STRATEGY_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>

#include <ANTL/XGCD/xgcd.hpp>
#include <ANTL/Arithmetic/mul_exact.hpp>

#define RED_CANTOR 0
#define RED_FAST 1

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class ReduceStrategy {
    protected:
      T Delta;
      T hx;
      long genus;
      bool is_init;

    public:
      ReduceStrategy() { is_init = false; };
      virtual ~ReduceStrategy() {};

      // Initialize field invariants
      void init(const T & Din, const T & hin, long gin=0) {
        if (is_init) {
          Delta.kill();
          hx.kill();
        }

        Delta = Din;
        hx = hin;
        genus = gin;
        is_init = true;
      }

      // Generic ideal reduction definition
      virtual void reduce(QuadraticIdealBase<T> & A) = 0;
  };

}// ANTL

#endif // guard
