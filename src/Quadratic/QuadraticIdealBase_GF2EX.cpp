/**
 * @file QuadraticIdealBase_GF2EX.cpp
 * @author Michael Jacobson
 * @remark Reduced quadratic ideal function specializations (GF2EX base type).
 */

#include <QuadraticIdealBase.hpp>

template <> 
void 
QuadraticIdealBase<GF2EX>::test_ideal (string msg)
{
  GF2EX temp = b * b + b*QO->getH() + a * c;

  if (temp != QO->getF())
    {
      cout << "ERROR " << msg << "!  wrong discriminant!" << endl;
      cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
      cout << "hx = " << QO->getH() << endl;
      cout << "fx = " << QO->getF() << endl;
      cout << "b^2 + bh + ac = " << temp << endl;
      exit (1);
    }
}



//
// QuadraticIdealBase<GF2EX>::assign_one()
//
// Task:
//      set to the unit ideal of the current quadratic_order
//

template <> 
void 
QuadraticIdealBase<GF2EX>::assign_one ()
{
  set (a);
  clear (b);
  c = QO->getF();
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
QuadraticIdealBase<GF2EX>::assign_prime (const GF2EX & p)
{
  GF2EX temp, Dp;
  long jac;

  if (!DetIrredTest (p))
    return false;

  jac = ressol (temp, QO->getH(), QO->getF(), p);

  if (jac < 0)
    return false;

  a = p;
  b = temp;

  // c = (b * b + b*hx +  fx) / a;
  sqr(c,b);
  mul(temp,b,QO->getH());
  add(c,c,temp);
  add(c,c,QO->getF());
  div(c,c,a);

  if (jac == 0)
    {
      // non-invertible if gcd(a,h,(f+bh+b^2)/a) <> 1
      temp = GCD(c,GCD(p,QO->getH()));
      if (!NTL::IsOne(temp))
	return false;
    }

  return true;
}



//
// conjugate()
//
// Task:
//      computes the conjugate of A
//

template <>
void
conjugate (QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A)
{
  GF2EX temp;

  C.a = A.a;
  C.b = A.b + QO->getH();
  if (deg(C.b) < deg(C.a))
    C.c = A.c;
  else {
    rem(C.b,C.b,C.a);

    // C.c = (C.b * C.b + C.b * hx +  fx) / C.a;
    sqr(C.c,C.b);
    mul(temp,C.b,QO->getH());
    add(C.c,C.c,temp);
    add(C.c,C.c,QO->getF());
    div(C.c,C.c,C.a);
  }
}
