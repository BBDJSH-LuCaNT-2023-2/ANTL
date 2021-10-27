/**
 * @file qo_multiply.hpp
 * @author Michael Jacobson
 * @brief Class defining abstract methods for ideal multiplication.
 */

#ifndef MULTIPLY_STRATEGY_H
#define MULTIPLY_STRATEGY_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>

//#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/XGCD/xgcd.hpp>
#include <ANTL/Arithmetic/mul_exact.hpp>

#define MUL_CANTOR 0
#define MUL_NUCOMP 1
#define MUL_EXPLICIT 2

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

  template <class T> class QuadraticIdealBase;

  template <class T> class MultiplyStrategy {
    protected:
      T    Delta;
      T    hx;
      long genus;
      bool is_init;

    public:
               MultiplyStrategy() {}
      virtual ~MultiplyStrategy() {}

      // Initialize field invariants
      void init(const T & Din, const T & hin, long gin=0) {
//         if (is_init) {
//           Delta.kill();
//           hx.kill();
//         }

        Delta = Din;
        hx = hin;
        genus = gin;
        is_init = true;
      }

      void getDelta() {
        std::cout << Delta << std::endl;
      }

      // Generic ideal multiplication definition
      virtual void multiply(QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> & A, const QuadraticIdealBase<T> & B) = 0;
  };

} // ANTL

#endif // guard
