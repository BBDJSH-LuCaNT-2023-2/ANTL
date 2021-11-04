/**
 * @file qo_square.hpp
 * @author Michael Jacobson
 * @brief Class defining abstract methods for ideal squaring.
 */

#ifndef SQUARE_STRATEGY_H
#define SQUARE_STRATEGY_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>

#include <ANTL/XGCD/xgcd.hpp>
#include <ANTL/Arithmetic/mul_exact.hpp>

#define SQR_CANTOR 0
#define SQR_NUCOMP 1
#define SQR_EXPLICIT 2

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class SquareStrategy {
    protected:
      T Delta;
      T hx;
      long genus;
      bool is_init;

    public:
      SquareStrategy() { is_init = false; };
      virtual ~SquareStrategy() {};

      // Initialize field invariants
      void init(const T & delta_in, const T & h_in, long g_in=0) {
//         if (is_init) {
//           Delta.kill();
//           hx.kill();
//         }

        Delta = delta_in;
        hx = h_in;
        genus = g_in;
        is_init = true;
      };

    // Generic ideal squaring definition
    virtual void square(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A) = 0;
  };

} //ANTL

#endif // guard
