/**
 * @file QuadraticIdealBase_ZZ.cpp
 * @author Michael Jacobson
 * @remark Primitive quadratic ideal function specializations (ZZ base type).
 */

#include <QuadraticIdealBase.hpp>


template <>
void
QuadraticIdealBase<ZZ>::test_ideal(string msg)
{
  ZZ tval = b*b - 4*a*c;
  if (tval != QO->getDiscriminant()) {
    cout << "ERROR " << msg << "!  wrong discriminant!" << endl;
    cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
    cout << "Delta = " << QO->getDiscriminant() << endl;
    cout << "b^2 - 4ac = " << tval << endl;
    exit(1);
  }
}



//
// QuadraticIdealBase<T>::assign_one()
//
// Task:
//      set to the unit ideal of the current quadratic_order
//

template <>
void
QuadraticIdealBase<ZZ>::assign_one()
{
  set(a);
  if (rem(QO->getDiscriminant(),4) == 1) {
    set(b);
    sub(c,1,QO->getDiscriminant());
  }
  else {
    clear(b);
    negate(c,QO->getDiscriminant());
  }
  div(c,c,a);
  RightShift(c,c,2);
}




//
// assign_prime()
//
// Task:
//      computes a reduced representative of the equivalence class containing
//      the ideal lying over the prime p.  If such an ideal doesnot exist,
//      false is returned.
//

template <>
bool
QuadraticIdealBase<ZZ>::assign_prime (const ZZ & p)
{
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
    else
      {
	long D8 = rem(QO->getDiscriminant(),8);

	if (D8 < 0)
	  add(D8,D8,8);

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
