/**
 * @file QuadraticIdealBase_long.cpp
 * @author Michael Jacobson
 * @remark Primitive quadratic ideal function specializations (long base type).
 */

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

using namespace ANTL;

template <> void QuadraticIdealBase<long>::ensure_valid(std::string msg) {
    long tval = b*b - 4*a*c;
    if (tval != QO->get_discriminant()) {
      cout << "ERROR " << msg << "!  wrong discriminant!" << endl;
      cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
      cout << "Delta = " << QO->get_discriminant() << endl;
      cout << "b^2 - 4ac = " << tval << endl;
      exit(1);
    }
  }

//
// QuadraticIdealBase<T>::assign_one()
//
// Task:
//      Set to the unit ideal of the current QuadraticOrder
//

template <> void QuadraticIdealBase<long>::assign_one() {
  a = 1;
  if ((QO->get_discriminant() & 3) == 1) {
    b = 1;
    c = 1 - QO->get_discriminant();
  }
  else {
    b = 0;
    c = -QO->get_discriminant();
  }
  c = c / (a << 2);
}

//
// assign_prime()
//
// Task:
//      computes a reduced representative of the equivalence class containing
//      the ideal lying over the prime p.  If such an ideal does not exist,
//      false is returned.
//

template <> bool QuadraticIdealBase<long>::assign_prime (const long & p) {

  long temp, Dp;
  long jac;

  if (!ProbPrime (p))
    return false;

  Dp = QO->get_discriminant() % p;
  if (Dp < 0)
    Dp += p;

  if (p == 2) {
    if (Dp == 0) {
      Dp = QO->get_discriminant() % 4;
      if (Dp < 0)  Dp += 4;
      if (Dp == 0) {
        Dp = (QO->get_discriminant() >> 2) % 4;
        if (Dp < 0)
          Dp += 4;
        if (Dp != 3)
          return false;
      }
      a = 2;
      b = 2;
      c = (4-QO->get_discriminant()) >> 3;
      return true;
    }
    else {
    long D8 = QO->get_discriminant() % 8;

    if (D8 < 0)
      D8 += 8;
    if (D8 == 1)
      jac = 1;
    else
      jac = -1;
    }
  }
  else {
    if (Dp == 0) {
      a = p;
      if (IsOdd (QO->get_discriminant()))
        b = p;
      else
        b = 0;
      c = (b*b - QO->get_discriminant()) / a;
      c >>= 2;

      temp = QO->get_discriminant() % (p * p);
      if (IsZero(to_ZZ(temp))) //An IsZero call was here, but there doesn't seem to be an implementation for that method yet...
    return false;
      else
        return true;
    }
    else
      jac = long(ANTL::Jacobi(to_ZZ(Dp), to_ZZ(p)));
  }

  if (jac < 0)
    return false;

  temp = to_long(SqrRootMod (to_ZZ(Dp), to_ZZ(p))); //SqrRootMod is not defined for input of type long.
  if (temp < 0)
    temp += p;
  if (IsOdd (QO->get_discriminant()) != IsOdd (temp))
    temp = p - temp;

  a = p;
  b = temp;
  c = (b * b - QO->get_discriminant()) / a;
  c >>= 2;

  return true;
}

// QuadraticIdealBase<T>::is_normal()
// Note: Not defined for positive definite forms.
// Task: tests if the ideal is normal.
template <> bool QuadraticIdealBase<long>::is_normal() {

  if(delta > 0) {
    if(abs(a) > floor_root_delta)
      return (-1*(abs(a)) < b && b <= abs(a));

    else
      return (floor_root_delta - 2*abs(a) < b && b <= floor_root_delta);
  }
}

// QuadraticIdealBase<T>::is_reduced()
//
// Task: tests if the ideal is reduced.
template <> bool QuadraticIdealBase<long>::is_reduced () {
  if(delta > 0) {
    long lbound, ubound, rootD;

    // lbound = abs(rootD - 2*abs(a))
    lbound = abs(floor_root_delta - 2*abs(a));

    // We assume the form is irrational, so there ought to be no case where delta is square
    ubound = floor_root_delta;

    return (lbound < b && b <= floor_root_delta);
  }

  // TODO: The case when delta < 0 remains untested!
  else if (delta < 0) {
    bool cond1 = ((abs(b) <= a) && (a <= c));
    bool cond2 = true;
    if( ((abs(b) == a) || (c == a)) && (b < 0))
      cond2 = false;

    return cond1 && cond2;
  };
}

template <> void QuadraticIdealBase<long>::normalize() {
static long a2, q, r, temp, nb;

  a2 = a << 1;
  if (a <= floor_root_delta) {
    temp = floor_root_delta - a2;
    if (!(temp < b && b <= floor_root_delta)) {
      temp = floor_root_delta - b;
      DivRem(q, r, temp, a2);

      nb = floor_root_delta - r;

      c += q * ((b + nb) >> 1);
      b = nb;
    }
  }
  else {
    if (!(b > -a && b <= a)) {
      DivRem(q, r, b, a2);

      if (r > a) {
        r -= a2;
        ++q;
      }

      c -= q * ((b + r) >> 1);
      b = r;
    }
  }
}
