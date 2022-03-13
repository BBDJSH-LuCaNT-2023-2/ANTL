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
  long delta, rootD;

  delta = b*b - 4*a*c;

  if(delta > 0) {
    // rootD = floor(sqrt(delta)) - [Recall  ANTL::SqrRoot(long a) = long floor(sqrt(a))]
    rootD = SqrRoot(abs(delta));

    if(abs(a) > rootD)
      return (-1*(abs(a)) < b && b <= abs(a));

    else
      return (rootD - 2*abs(a) < b && b <= rootD);
  }
}

// QuadraticIdealBase<T>::is_reduced()
//
// Task: tests if the ideal is reduced.
template <> bool QuadraticIdealBase<long>::is_reduced () {
  long delta;

  delta = b*b - 4*a*c;

  if(delta > 0) {
    long lbound, ubound, rootD;

    // rootD = floor(sqrt(delta)) - [Recall  ANTL::SqrRoot(long a) = long floor(sqrt(a))]
    rootD = SqrRoot(abs(delta));

    // lbound = abs(rootD - 2*abs(a))
    lbound = abs(rootD - 2*abs(a));

    // We assume the form is irrational, so there ought to be no case where delta is square
    ubound = rootD;

    return (lbound < b && b <= rootD);
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
  static long a2, delta, rootDelta, s;

  delta = b*b - 4*a*c;

  rootDelta = SqrRoot(delta);

  if(a <= rootDelta) {
    a2 = 2*abs(a);

    // Computing s, the normalizing integer,  per [BV07, pg. 108]
    s = sign(a)*((rootDelta - b)/a2);

    //c = a*s^2 + b*s + c
    c = a*s*s + b*s + c;

    //b = b + 2sa
    b = b + 2*s*a;
  }

  else {
    a2 = 2*abs(a);

    // Computing s, the normalizing integer,  per [BV07, pg. 108]
    s = sign(a)*((abs(a) - b)/a2);

    //c = a*s^2 + b*s + c
    c = a*s*s + b*s + c;

    //b = b + 2sa
    b = b + 2*s*a;
  }
}
