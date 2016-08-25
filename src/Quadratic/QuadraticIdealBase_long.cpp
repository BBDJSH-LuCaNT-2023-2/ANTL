/**
 * @file QuadraticIdealBase_long.cpp
 * @author Michael Jacobson
 * @remark Primitive quadratic ideal function specializations (long base type).
 */

#include <QuadraticIdealBase.hpp>


template <>
void
QuadraticIdealBase<long>::test_ideal(string msg)
{
  long tval = b*b - 4*a*c;
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
QuadraticIdealBase<long>::assign_one()
{
  a = 1;
  if ((QO->getDiscriminant() & 3) == 1) {
    b = 1;
    c = 1 - QO->getDiscriminant();
  }
  else {
    b = 0;
    c = -QO->getDiscriminant();
  }
  c = c / (a << 2);
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
QuadraticIdealBase<long>::assign_prime (const long & p)
{
  long temp, Dp;
  long jac;

  if (!ProbPrime (p))
    return false;

  Dp = QO->getDiscriminant() % p;
  if (Dp < 0)
    Dp += p;

  if (p == 2) {
    if (Dp == 0) {
      Dp = QO->getDiscriminant() % 4;
      if (Dp < 0)  Dp += 4;
      if (Dp == 0) {
        Dp = (QO->getDiscriminant() >> 2) % 4;
        if (Dp < 0)
          Dp += 4;
        if (Dp != 3)
          return false;
      }
      a = 2;
      b = 2;
      c = (4-QO->getDiscriminant()) >> 3;
      return true;
    }
    else
      {
	long D8 = QO->getDiscriminant() % 8;

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
      if (IsOdd (QO->getDiscriminant()))
        b = p;
      else
        b = 0;
      c = (b*b - QO->getDiscriminant()) / a;
      c >>= 2;

      temp = QO->getDiscriminant() % (p * p);
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
    temp += p;
  if (IsOdd (QO->getDiscriminant()) != IsOdd (temp))
    temp = p - temp;

  a = p;
  b = temp;
  c = (b * b - QO->getDiscriminant()) / a;
  c >>= 2;

  return true;
}
