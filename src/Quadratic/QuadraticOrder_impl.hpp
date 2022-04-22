/**
 * @file QuadraticOrder_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from QuadraticOrder.hpp only.
 */

#include <iostream>

namespace ANTL {
//
// constructor
//     - sets resets all the class members which are vectors
//

template <class T> QuadraticOrder<T>::QuadraticOrder(const T &D) {
  // test whether D is a valid discriminant
  if (deg(D) > !3) {
    // assign new values
    Delta = D;
    clear(hx);

    // compute genus
    g = deg(Delta);
    if (g & 1)
      g = (g - 1) >> 1;
    else {
      g = (g >> 1) - 1;
      //        if (qi_pair<T>::get_rootD()*qi_pair<T>::get_rootD() == Delta)
      // return false;
    }

    /*
    // initialize invariant values
    clear (R);
    clear (h);
    CL.clear();
    gens.clear();
    clear (L);
    */
    //Lfunc.init(Delta, QUADRATIC_MODE);
  }
}

//
// destructor
//

template <class T> QuadraticOrder<T>::~QuadraticOrder() {}

/*
//
// QuadraticOrder<T>::verbose()
//
// Task:
//      sets the verbosity of commands.  Currently, the following levels are
//      supported:
//         0 - nothing
//         1 - run times and some run-time data (data only for subexp algs.)
//         >1 - state information for subexponential algorithms
//

template < class T >
void QuadraticOrder<T>::verbose (long state)
{
  if (state <= 0)
    {
      info = 0;
    }
  else
    {
      info = state;
      L_function<T>::verbose(state);
    }
}
*/

// Multiply Strategy Setters
template <class T> void QuadraticOrder<T>::set_red_best(ReduceStrategy<T> &A) {
  A.init(Delta, hx, g);
  red_best = &A;
}

template <class T>
void QuadraticOrder<T>::set_red_plain_imag(ReducePlainImag<T> &A) {
  A.init(Delta, hx, g);
  red_plain_imag = &A;

  // If red_best is not set, or if Plain is the preferred strategy, set red_best
  // to A as well
  if (red_best == nullptr || (preferred_red == 0 && this->is_imaginary())) {
    red_best = &A;
  }
}

template <class T>
void QuadraticOrder<T>::set_red_plain_real(ReducePlainReal<T> &A) {
  A.init(Delta, hx, g);
  red_plain_real = &A;

  // If red_best is not set, or if Plain is the preferred strategy, set red_best
  // to A as well
  if (red_best == nullptr || (preferred_red == 0 && this->is_real())) {
    red_best = &A;
  }
}

template <class T> void QuadraticOrder<T>::set_red_fast(ReduceFast<T> &A) {
  A.init(Delta, hx, g);
  red_fast = &A;

  // If red_best is not set, or if Fast is the preferred strategy, set red_best
  // to A as well
  if (red_best == nullptr || preferred_red == 1) {
    red_best = &A;
  }
}

// Multiply Strategy Setters
template <class T>
void QuadraticOrder<T>::set_mul_best(MultiplyStrategy<T> &A) {
  mul_best = &A;
}

template <class T>
void QuadraticOrder<T>::set_mul_comp(MultiplyComp<T> &A) {
  mul_comp = &A;

  // If mul_best is not set, or if Nucomp is the preferred strategy, set
  // mul_best to A as well
  if (mul_best == nullptr || preferred_mul == 2) {
    mul_best = &A;
  }
}


template <class T> void QuadraticOrder<T>::set_mul_plain(MultiplyPlain<T> &A) {
  mul_plain = &A;

  // If mul_best is not set, or if Plain is the preferred strategy, set mul_best
  // to A as well
  if (mul_best == nullptr || preferred_mul == 0) {
    mul_best = &A;
  }
}

template <class T>
void QuadraticOrder<T>::set_mul_nucomp(MultiplyNucomp<T> &A) {
  mul_nucomp = &A;

  // If mul_best is not set, or if Nucomp is the preferred strategy, set
  // mul_best to A as well
  if (mul_best == nullptr || preferred_mul == 1) {
    mul_best = &A;
  }
}

// Square Strategy Setters
template <class T> void QuadraticOrder<T>::set_sqr_best(SquareStrategy<T> &A) {
  sqr_best = &A;
}

template <class T> void QuadraticOrder<T>::set_sqr_plain(SquarePlain<T> &A) {
  sqr_plain = &A;

  // If sqr_best is not set, or if Plain is the preferred strategy, set sqr_best
  // to A as well
  if (sqr_best == nullptr || preferred_sqr == 0) {
    sqr_best = &A;
  }
}

template <class T> void QuadraticOrder<T>::set_sqr_nudupl(SquareNudupl<T> &A) {
  sqr_nudupl = &A;

  // If sqr_best is not set, or if Nudupl is the preferred strategy, set
  // sqr_best to A as well
  if (sqr_best == nullptr || preferred_sqr == 0) {
    sqr_best = &A;
  }
}

// Cube Strategy Setters
template <class T> void QuadraticOrder<T>::set_cube_best(CubeStrategy<T> &A) {
  cube_best = &A;
}

template <class T> void QuadraticOrder<T>::set_cube_plain(CubePlain<T> &A) {
  cube_plain = &A;

  // If cube_best is not set, or if Plain is the preferred strategy, set
  // cube_best to A as well
  if (cube_best == nullptr || preferred_cube == 0) {
    cube_best = &A;
  }
}

template <class T> void QuadraticOrder<T>::set_cube_nucube(CubeNucube<T> &A) {
  cube_nucube = &A;

  // If cube_best is not set, or if Nucube is the preferred strategy, set
  // cube_best to A as well
  if (cube_best == nullptr || preferred_cube == 1) {
    cube_best = &A;
  }
}

template <class T> void QuadraticOrder<T>::set_cube_mulsqr(CubeMulSqr<T> &A) {
  cube_mulsqr = &A;

  // If cube_best is not set, or if Mulsqr is the preferred strategy, set
  // cube_best to A as well
  if (cube_best == nullptr || preferred_cube == 2) {
    cube_best = &A;
  }
}

template <class T> ReduceStrategy<T> *QuadraticOrder<T>::get_red_best() {
  return red_best;
}
template <class T> ReducePlainImag<T> *QuadraticOrder<T>::get_red_plain_imag() {
  return red_plain_imag;
}
template <class T> ReducePlainReal<T> *QuadraticOrder<T>::get_red_plain_real() {
  return red_plain_real;
}
template <class T> ReduceFast<T> *QuadraticOrder<T>::get_red_fast() {
  return red_fast;
}

template <class T> MultiplyStrategy<T> *QuadraticOrder<T>::get_mul_best() {
  return mul_best;
}
template <class T> MultiplyComp<T> *QuadraticOrder<T>::get_mul_comp() {
  return mul_comp;
}
template <class T> MultiplyPlain<T> *QuadraticOrder<T>::get_mul_plain() {
  return mul_plain;
}
template <class T> MultiplyNucomp<T> *QuadraticOrder<T>::get_mul_nucomp() {
  return mul_nucomp;
}

template <class T> SquareStrategy<T> *QuadraticOrder<T>::get_sqr_best() {
  return sqr_best;
}
template <class T> SquarePlain<T> *QuadraticOrder<T>::get_sqr_plain() {
  return sqr_plain;
}
template <class T> SquareNudupl<T> *QuadraticOrder<T>::get_sqr_nudupl() {
  return sqr_nudupl;
}

template <class T> CubeStrategy<T> *QuadraticOrder<T>::get_cube_best() {
  return cube_best;
}
template <class T> CubePlain<T> *QuadraticOrder<T>::get_cube_plain() {
  return cube_plain;
}
template <class T> CubeNucube<T> *QuadraticOrder<T>::get_cube_nucube() {
  return cube_nucube;
}
template <class T> CubeMulSqr<T> *QuadraticOrder<T>::get_cube_mulsqr() {
  return cube_mulsqr;
}

//
// QuadraticOrder<T>::IsEqual()
//
// Task:
//      tests if the function fields are equal
//

template <class T>
bool QuadraticOrder<T>::IsEqual(const QuadraticOrder<T> &QO) const {
  return (Delta == QO.Delta && hx == QO.hx);
}

//
// operator ==
//
// Task:
//      tests if QO1 and QO2 are equal (same discriminant)
//

template <class T>
bool operator==(const QuadraticOrder<T> &QO1, const QuadraticOrder<T> &QO2) {
  return (QO1.IsEqual(QO2));
}

//
// operator !=
//
// Task:
//      tests if QO1 and QO2 are not equal
//

template <class T>
bool operator!=(const QuadraticOrder<T> &QO1, const QuadraticOrder<T> &QO2) {
  return (!QO1.IsEqual(QO2));
}

//
// QuadraticOrder<T>::is_imaginary()
//
// Task:
//      returns true if the function field is imaginary
//

template <class T> bool QuadraticOrder<T>::is_imaginary() const {
  return (deg(Delta) & 1);
}

//
// QuadraticOrder<T>::IsUnusual()
//
// Task:
//      returns true if the function field is real
//

template <class T> bool QuadraticOrder<T>::IsUnusual() const {
  return (!(deg(Delta) & 1) && !test_Dcoeff(Delta));
}

//
// QuadraticOrder<T>::is_real()
//
// Task:
//      returns true if the function field is real
//

template <class T> bool QuadraticOrder<T>::is_real() const {
  return (!(deg(Delta) & 1) && test_Dcoeff(Delta));
}

/*
//
// randomImaginaryOrder
//
// Task:
//      generates a random imaginary order
//

template < class T >
QuadraticOrder<T> & randomImaginaryOrder (long size, bool prime)
{
  // select a random discriminant with given genus
  T disc;

  do
    {
      random (disc, 2 * size + 2);
      MakeMonic (disc);
    }
  while (deg (disc) != (2 * size + 1) || !DetIrredTest (disc));

  return new QuadraticOrder<T>(disc);
}



//
// randomUnusualOrder
//
// Task:
//      generates a random imaginary order
//

template < class T >
QuadraticOrder<T> &
randomUnusualOrder (long size, bool prime)
{
  // select a random discriminant with given genus

  T disc;

  do
    {
      random (disc, 2 * size + 3);
      MakeMonic (disc);
    }
  while (deg (disc) != (2 * size + 2) || test_Dcoeff(disc));

  return new QuadraticOrder<T>(disc);
}



//
// randomRealOrder
//
// Task:
//      generates a random real order
//

template < class T >
QuadraticOrder<T> &
randomRealOrder (long size, bool prime)
{
  // select a random discriminant with given genus

  T disc;

  do
    {
      random (disc, 2 * size + 3);
      MakeMonic (disc);
    }
  while (deg (disc) != (2 * size + 2) || !DetIrredTest (disc));

  return new QuadraticOrder<T>(disc);
}
*/

/*
//
// read_from_file
//
// Task:
//      inputs a QuadraticOrder from the istream in.
//

template < class T >
void
QuadraticOrder<T>::read_from_file(std::istream & in)
{
  T newD;

  in >> newD;
  assign (newD);
}
*/

//
// operator >>
//
// Task:
//      inputs a QuadraticOrder from the istream in.
//

template <class T>
std::istream &operator>>(std::istream &in, QuadraticOrder<T> &QO) {
  QO.read_from_file(in);
  return in;
}

//
// operator <<
//
// Task:
//      outputs a QuadraticOrder to the ostream out.
//

/*
template < class T >
void QuadraticOrder<T>::write_to_file(std::ostream & out) const
{
  out << "Quadratic order (genus " << g << "):" << endl;
  output_FF (out, Delta);

  if (!is_zero ())
    {
      if (is_R_computed () && is_real ())
        out << "   R = " << R << endl;
      if (is_h_computed ())
        out << "   h = " << h << endl;
      if (is_CL_computed ())
        out << "   CL = " << CL << endl;
      if (rank > 0)
        {
          out << "   generators:" << endl;
          for (long i = 0; i < rank; ++i)
            {
              out << "      " << gens[i];
              if (i < rank - 1)
                out << ",";
              out << endl;
            }
        }

      if (is_L_computed ())
        out << "   L(1) = " << L << endl;
    }
}
*/

//
//
// operator <<
//
// Task:
//      outputs a QuadraticOrder to the ostream out.
//

template <class T>
std::ostream &operator<<(std::ostream &out, const QuadraticOrder<T> &QO) {
  QO.write_to_file(out);
  return out;
}

} // namespace ANTL
