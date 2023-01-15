#ifndef MULTIPLY_NUCOMP_OPT_H
#define MULTIPLY_NUCOMP_OPT_H

#include <ANTL/Quadratic/Multiply/MultiplyStrategy.hpp>
#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

NTL_CLIENT
using namespace ANTL;

namespace ANTL {

template <class T> class MultiplyStrategy;
template <class T> class QuadraticIdealBase;
template <class T> class QuadraticNumber;

template <class T> class MultiplyNucompOpt : public MultiplyStrategy<T> {

  using MultiplyStrategy<T>::Delta;
  using MultiplyStrategy<T>::hx;
  using MultiplyStrategy<T>::genus;
  using MultiplyStrategy<T>::is_init;

protected:
  QuadraticNumber<T> *RelativeGenerator;
  ZZ NC_BOUND; // termination bound for NUCOMP = floor(|D|^1/4)

public:

  MultiplyNucompOpt(QuadraticOrder<T> &inQO);

  ~MultiplyNucompOpt(){};

  void init(const T &delta_in, const T &h_in, long g_in = 0) {
    MultiplyStrategy<T>::init(delta_in, h_in, g_in);
  };

  QuadraticNumber<T> *get_RelativeGenerator() { return RelativeGenerator; }

  void set_RelativeGenerator(QuadraticNumber<T> &QN) {
    RelativeGenerator = &QN;
  }

  //     nucomp();
  //     Task:
  //          computes a reduced ideal equivalent to the product of two ideals
  //          using the NUCOMP algorithm of Shanks.

  void multiply(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A,
                const QuadraticIdealBase<T> &B);

private:
  void construct_relative_generator(T &rel_gen_a, T &rel_gen_b, T &rel_gen_d,
                                    QuadraticIdealBase<T> &C, T OB, T BB,
                                    T S);
};

//
// Declare specialized methods
//

template <>
void MultiplyNucompOpt<ZZ>::init(const ZZ &delta_in, const ZZ &h_in, long g_in);
template <>
void MultiplyNucompOpt<long>::init(const long &delta_in, const long &h_in,
                                   long g_in);

template <>
void MultiplyNucompOpt<ZZ>::multiply(QuadraticIdealBase<ZZ> &C,
                                     const QuadraticIdealBase<ZZ> &A,
                                     const QuadraticIdealBase<ZZ> &B);

template <>
void MultiplyNucompOpt<long>::multiply(QuadraticIdealBase<long> &C,
                                       const QuadraticIdealBase<long> &A,
                                       const QuadraticIdealBase<long> &B);

template <>
void MultiplyNucompOpt<GF2EX>::multiply(QuadraticIdealBase<GF2EX> &C,
                                        const QuadraticIdealBase<GF2EX> &A,
                                        const QuadraticIdealBase<GF2EX> &B);

} // namespace ANTL

// Unspecialized template definitions.
#include "../../../../src/Quadratic/Multiply/MultiplyNucomp_Opt_impl.hpp"

#endif // guard
