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

  if (rem(QO->getDiscriminant(),4) == 1) {
    set(b);
    sub(c,1,QO->getDiscriminant());
  }

  else {
    clear(b);
    NTL::negate(c,QO->getDiscriminant());
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

  rem(Dp,QO->getDiscriminant(),p);
  if (Dp < 0)
    add(Dp,Dp,p);

  if (p == 2) {
    if (IsZero(Dp)) {
      rem(QO->getDiscriminant(),4);
        if (Dp < 0)
          add(Dp,Dp,4);

      if (IsZero(Dp)) {
        RightShift(temp,QO->getDiscriminant(),2);
        Dp = rem(temp,4);
        if (Dp < 0)
          add(Dp,Dp,4);
        if (Dp != 3)
          return false;
      }
      a = 2;
      b = 2;
      sub(c,4,QO->getDiscriminant());
      RightShift(c,c,3);
      return true;
    }

    else {
      ZZ D8;
      rem(D8, QO->getDiscriminant(),ZZ(8));

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
      if (IsOdd (QO->getDiscriminant()))
        b = p;
      else
        b = 0;
      sqr(c,b);
      sub(c,c,QO->getDiscriminant());
      div(c,c,a);
      RightShift(c,c,2);

      rem(temp,QO->getDiscriminant(),p*p);
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

  if (IsOdd (QO->getDiscriminant()) != IsOdd (temp))
    sub(temp,p,temp);

  a = p;
  b = temp;
  sqr(c,b);
  sub(c,c,QO->getDiscriminant());
  div(c,c,a);
  RightShift(c,c,2);

  return true;
}

// QuadraticIdealBase<T>::is_reduced()
//
// Task: tests if the ideal is reduced.
template <> bool QuadraticIdealBase<ZZ>::is_reduced () {
  ZZ D_ZZ = b^2 - 4*a*c;
  long D = to_long(D_ZZ);

  if(D < 0) {
    bool cond1 = ((abs(b) <= a) && (a <= c));
    bool cond2 = true;
    if( ((abs(b) == a) || (c == a)) && (b < 0))
      cond2 = false;

    return cond1 && cond2;
  }

  else if (D > 0) {
    long lbound = abs(std::sqrt(D) - to_long(2*abs(c)));
    long ubound = std::sqrt(D);

    return (lbound < b && b < ubound);
  };
}
