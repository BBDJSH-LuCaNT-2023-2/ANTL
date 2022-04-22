/**
 * @file QuadraticIdealBase_ZZ.cpp
 * @author Michael Jacobson
 * @remark Primitive quadratic ideal function specializations (ZZ base type).
 */

#include <ANTL/Quadratic/QuadraticIdealBase.hpp>

using namespace ANTL;

// QuadraticIdealBase<T>::assign_one()
//
// Task: set to the unit ideal of the current quadratic_order

template <> void QuadraticIdealBase<ZZ>::assign_one() {
  set(a);

  if (rem(QO->get_discriminant(),4) == 1) {
    set(b);
    sub(c,1,QO->get_discriminant());
  }

  else {
    clear(b);
    NTL::negate(c,QO->get_discriminant());
  }

  div(c,c,a);
  RightShift(c,c,2);
}

// assign_prime()
//
// Task: Computes an ideal lying over the prime p. If such an ideal does not exist, false is returned.

template <> bool QuadraticIdealBase<ZZ>::assign_prime (const ZZ & p) {
  ZZ temp, Dp;
  long jac;

  if (!ProbPrime (p))
    return false;

  rem(Dp,QO->get_discriminant(),p);
  if (Dp < 0)
    add(Dp,Dp,p);

  if (p == 2) {
    if (IsZero(Dp)) {
      rem(QO->get_discriminant(),4);
        if (Dp < 0)
          add(Dp,Dp,4);

      if (IsZero(Dp)) {
        RightShift(temp,QO->get_discriminant(),2);
        Dp = rem(temp,4);
        if (Dp < 0)
          add(Dp,Dp,4);
        if (Dp != 3)
          return false;
      }
      a = 2;
      b = 2;
      sub(c,4,QO->get_discriminant());
      RightShift(c,c,3);
      return true;
    }

    else {
      ZZ D8;
      rem(D8, QO->get_discriminant(),ZZ(8));

      if (D8 < 0)
        add(D8,D8,ZZ(8));

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
      sqr(c,b);
      sub(c,c,QO->get_discriminant());
      div(c,c,a);
      RightShift(c,c,2);

      rem(temp,QO->get_discriminant(),p*p);
      if (IsZero (temp))
        return false;
      else
        return true;
    }
    else
      jac = Jacobi_base(Dp,p);
  }

  if (jac < 0)
    return false;

  temp = SqrRootMod (Dp, p);
  if (temp < 0)
    add(temp,temp,p);

  if (IsOdd (QO->get_discriminant()) != IsOdd (temp))
    sub(temp,p,temp);

  a = p;
  b = temp;
  sqr(c,b);
  sub(c,c,QO->get_discriminant());
  div(c,c,a);
  RightShift(c,c,2);

  return true;
}

// QuadraticIdealBase<T>::is_normal()
// Note: Not defined for positive definite forms.
// Task: tests if the ideal is normal.
template <> bool QuadraticIdealBase<ZZ>::is_normal() {
  ZZ Delta, FloorRootDelta;

  if (a < 0) {
    return false;
  }

  Delta = this->get_QO()->get_discriminant();
  FloorRootDelta = SqrRoot(Delta);

  if(FloorRootDelta < a) {
    return (-a < b) && (b <= a);
  }

  else {
    return (FloorRootDelta - 2*a < b) && (b < FloorRootDelta + 1);
  }
}

// QuadraticIdealBase<T>::is_reduced()
//
// Task: tests if the ideal is reduced.
template <> bool QuadraticIdealBase<ZZ>::is_reduced() {
  ZZ Delta, FloorRootDelta;

  Delta = this->get_QO()->get_discriminant();

  if(Delta > 0) {
    FloorRootDelta = SqrRoot(Delta);

    // Below is a check for |sqrt(delta) - 2|a|| < b < sqrt(delta)
    // It is slightly modified, so it can use floor(sqrt(delta)) instead
    return (abs(FloorRootDelta - 2*a) < b) && (b < FloorRootDelta + 1);
  }

  // TODO: The case when delta < 0 remains untested!
  else if (Delta < 0) {
    bool cond1 = ((abs(b) <= a) && (a <= c));
    bool cond2 = true;

    if( ((abs(b) == a) || (c == a)) && (b < 0)) {
      cond2 = false;
    }

    return cond1 && cond2;
  };

}

template <> void QuadraticIdealBase<ZZ>::normalize() {
  ZZ Delta, FloorRootDelta, q, r;

  Delta = this->get_QO()->get_discriminant();
  FloorRootDelta = SqrRoot(Delta);

  abs(a, a);

  if(FloorRootDelta < a) {
    b += ((a-b) / (2*a)) * (2*a);
  }

  else {
    b += ((FloorRootDelta-b) / (2*a)) * (2*a);
  }

  c = (Delta - (b*b))/(-4*a);

}
