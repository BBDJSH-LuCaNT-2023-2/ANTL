#ifndef SQUARE_NUDUPL_OPT_H
#define SQUARE_NUDUPL_OPT_H

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>
#include <ANTL/Quadratic/Square/SquareStrategy.hpp>

NTL_CLIENT

using namespace ANTL;

namespace ANTL {

template <class T> class SquareStrategy;
template <class T> class QuadraticIdealBase;
template <class T> class QuadraticNumber;

template <class T> class SquareNuduplOpt : public SquareStrategy<T> {
  using SquareStrategy<T>::Delta;
  using SquareStrategy<T>::hx;
  using SquareStrategy<T>::genus;
  using SquareStrategy<T>::is_init;

protected:
  QuadraticNumber<T> *RelativeGenerator;
  ZZ NC_BOUND; // termination bound for NUCOMP = floor(|D|^1/4)

public:
  ~SquareNuduplOpt(){};

  void init(const T &delta_in, const T &h_in, long g_in = 0) {
    SquareStrategy<T>::init(delta_in, h_in, g_in);
  };

  QuadraticNumber<T> *get_RelativeGenerator() { return RelativeGenerator; }

  void set_RelativeGenerator(QuadraticNumber<T> &QN) {
    RelativeGenerator = &QN;
  }

  // nudupl
  void square(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);

private:
  void construct_relative_generator(T &rel_gen_a, T &rel_gen_b, T &rel_gen_d,
                                    QuadraticIdealBase<T> &C, T OB, T BB,
                                    T S);
};

// Declare specialized methods
template <>
void SquareNuduplOpt<ZZ>::init(const ZZ &delta_in, const ZZ &h_in, long g_in);
template <>
void SquareNuduplOpt<ZZ>::square(QuadraticIdealBase<ZZ> &C,
                              const QuadraticIdealBase<ZZ> &A);

template <>
void SquareNuduplOpt<long>::init(const long &delta_in, const long &h_in,
                              long g_in);
template <>
void SquareNuduplOpt<long>::square(QuadraticIdealBase<long> &C,
                                const QuadraticIdealBase<long> &A);

template <>
void SquareNuduplOpt<GF2EX>::square(QuadraticIdealBase<GF2EX> &C,
                                 const QuadraticIdealBase<GF2EX> &A);

} // namespace ANTL

// Unspecialized template definitions.

// TODO: impl file
// #include "../src/Quadratic/Square/SquareNuduplOpt_impl.hpp"

#endif // guard
