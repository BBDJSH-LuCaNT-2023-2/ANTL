/**
 * @file L_function_util.cpp
 * @author Leonard Nooy
 * @version $Header$
 */

#include <ANTL/common.hpp>
#include <ANTL/L_function/L_function_util.hpp>

namespace ANTL
{

  void Bach_Table (long Q, double &A, double &B)
  {
    if (Q >= 1000000)
      {
	A = 6.246;
	B = 16.217;
      }
    else if (Q >= 500000)
      {
	A = 6.269;
	B = 16.409;
      }
    else if (Q >= 100000)
      {
	A = 6.338;
	B = 17.031;
      }
    else if (Q >= 50000)
      {
	A = 6.378;
	B = 17.397;
      }
    else if (Q >= 10000)
      {
	A = 6.510;
	B = 18.606;
      }
    else if (Q >= 5000)
      {
	A = 6.593;
	B = 19.321;
      }
    else if (Q >= 1000)
      {
	A = 6.897;
	B = 21.528;
      }
    else if (Q >= 500)
      {
	A = 7.106;
	B = 22.845;
      }
    else if (Q >= 100)
      {
	A = 7.962;
	B = 27.145;
      }
    else if (Q >= 50)
      {
	A = 8.628;
	B = 29.587;
      }
    else if (Q >= 10)
      {
	A = 12.170;
	B = 38.831;
      }
    else
      {
	A = 16.397;
	B = 47.183;
      }
  }

  bool isInteger (CC < RR > &X, double diff)
  {
    RR Tr, Ti, temp1, temp2;

    Tr = X.real ();
    Ti = X.imaginary ();
    temp1 = abs (Tr - round (Tr));
    temp2 = abs (Ti - round (Ti));
    return (temp1 < diff && temp2 < diff);
  }

  void testInteger (CC < RR > &X, RR & diff)
  {
    RR Tr, Ti, temp1, temp2;

    Tr = X.real ();
    Ti = X.imaginary ();
    temp1 = abs (Tr - round (Tr));
    temp2 = abs (Ti - round (Ti));
    diff = sqrt(temp1*temp1 + temp2*temp2);
#ifdef EXTRA
    cout << "TEST:  " << X << ", t1 = " << temp1 << ", t2 = " << temp2 << ", d = " << diff << endl;
#endif

    /*
      if (temp1 > temp2)
      diff = temp1;
      else
      diff = temp2;
    */
  }

  void Cornacchia (ZZ & p, ZZ & x, ZZ & y)
  {
    ZZ x0, a, b, l, r, p4;

    p4 = p << 2;

    SqrRootMod (x0, p - 4, p);
    if (IsOdd (x0))
      x0 = p - x0;

    a = p << 1;
    b = x0;
    SqrRoot (l, p4);

    while (b > l)
      {
	rem (r, a, b);
	a = b;
	b = r;
      }

    sqr (l, b);
    r = (p4 - l) >> 2;
    x = b >> 1;
    SqrRoot (y, r);

    if (IsOdd (y))
      swap (x, y);

    long m8 = rem (x, 4);

    if (m8 == 1)
      x = -x;

    m8 = rem (y, 4);
    if (m8 != 2)
      y = -y;

    m8 = rem (x * y, 8);
    if (m8 < 0)
      m8 += 8;
    if (m8 != 2)
      y = -y;
  }

  void Cornacchia (long &p, long &x, long &y)
  {
    ZZ XX, p4;
    long x0, a, b, l, r;

    p4 = to_ZZ (p) << 2;

    SqrRootMod (XX, to_ZZ (p - 4), to_ZZ (p));
    conv (x0, XX);
    if (x0 & 1)
      x0 = p - x0;

    a = p << 1;
    b = x0;
    conv (l, SqrRoot (p4));

    while (b > l)
      {
	r = a % b;
	a = b;
	b = r;
      }

    sqr (XX, to_ZZ (b));
    conv (r, (p4 - XX) >> 2);
    x = b >> 1;
    y = SqrRoot (r);

    if (y & 1)
      {
	r = x;
	x = y;
	y = r;
      }

    if ((x & 3) == 1)
      x = -x;

    if ((y & 3) != 2)
      y = -y;

    long m8 = (x * y) & 7;

    if (m8 < 0)
      m8 += 8;
    if (m8 != 2)
      y = -y;
  }

//   void Cornacchia (long long &p, long long &x, long long &y)
//   {
//     ZZ XX, p4;
//     long long x0, a, b, l, r;
//
//     p4 = to<ZZ> (p) << 2;
//
//     SqrRootMod (XX, to<ZZ> (p - 4), to<ZZ> (p));
//     x0 = to<long long>(XX);
//     if (x0 & 1)
//       x0 = p - x0;
//
//     a = p << 1;
//     b = x0;
//     l = to<long long>(SqrRoot (p4));
//
//     while (b > l)
//       {
// 	r = a % b;
// 	a = b;
// 	b = r;
//       }
//
//     sqr (XX, to<ZZ> (b));
//     r = to<long long>((p4 - XX) >> 2);
//     x = b >> 1;
//     y = SqrRoot (r);
//
//     if (y & 1)
//       {
// 	r = x;
// 	x = y;
// 	y = r;
//       }
//
//     if ((x & 3) == 1)
//       x = -x;
//
//     if ((y & 3) != 2)
//       y = -y;
//
//     long m8 = (x * y) & 7;
//
//     if (m8 < 0)
//       m8 += 8;
//     if (m8 != 2)
//       y = -y;
//   }


  int
  compute_r2(const ZZ_pX & D)
  {
    long j;
    vec_pair_ZZ_pX_long facts;

    CanZass(facts,D);
    int r2 = 0;
    for (j=0; j<facts.length(); ++j)
      if (deg(facts[j].a) == 1) ++r2;
    return r2;
  }

  int
  compute_r2(const zz_pX & D)
  {
    long j;
    vec_pair_zz_pX_long facts;

    CanZass(facts,D);
    int r2 = 0;
    for (j=0; j<facts.length(); ++j)
      if (deg(facts[j].a) == 1) ++r2;
    return r2;
  }

  int
  compute_r2(const ZZ_pEX & D)
  {
    long j;
    vec_pair_ZZ_pEX_long facts;

    CanZass(facts,D);
    int r2 = 0;
    for (j=0; j<facts.length(); ++j)
      if (deg(facts[j].a) == 1) ++r2;
    return r2;
  }

  int
  compute_r2(const zz_pEX & D)
  {
    long j;
    vec_pair_zz_pEX_long facts;

    CanZass(facts,D);
    int r2 = 0;
    for (j=0; j<facts.length(); ++j)
      if (deg(facts[j].a) == 1) ++r2;
    return r2;
  }

  int
  compute_r2(const GF2EX & D)
  {
    long j;
    vec_pair_GF2EX_long facts;
    //MDV In case D's not monic..
    GF2EX D2 = D;
    MakeMonic(D2);
	
    CanZass(facts,D2);
    int r2 = 0;
    for (j=0; j<facts.length(); ++j)
      if (deg(facts[j].a) == 1) ++r2;
    return r2;
  }

  double
  eval_c(const RR & D) {
    double c;

    c = double(2) * ( double(5)/double(3) + to_double(log(D))) / double(3);
    return c;
  }


  double
  eval_U(double Q) {
    double temp,U;

    temp = double(2.0)*Q - double(1);
    U = double(0.5)*(temp*temp*(std::log(temp) - double(0.5)) - Q*Q*(std::log(Q) - double(0.5)));

    return U;
  }


  double
  eval_G(double Q, double lQ, double U, double l) {
    double G;

    G = double(1) + std::sqrt(double(8)) + double(5)*l + double(7)*l / lQ + double(4)*l / (lQ*lQ);
    G = std::sqrt(Q*Q*Q)*G / U;

    G += double(12) / (Q*lQ);
    G += double(8) / (Q*lQ*lQ);
    G += double(4) / (Q*lQ*lQ*lQ);
    G += double(12) / (Q*Q*lQ);
    G += double(15) / (double(2)*Q*Q*lQ*lQ);
    G += double(3) / (Q*Q*lQ*lQ*lQ);

    return G;
  }


  double
  eval_H(double Q, double lQ, double U, double l) {
    double H;
    double rQ = std::sqrt(Q*Q*Q);

    H = double(4)*double(1.25506)*l*rQ / U;
    H += double(3)*double(1.25506) / (std::log(double(2))*rQ);
    H += double(6) / Q;
    H += (double(10) + std::log(double(4))) / (Q*lQ);
    H += double(6) / (Q*lQ*lQ);
    H += double(2) / (Q*lQ*lQ*lQ);
    H += double(4) / (Q*Q*lQ);
    H += double(5) / (double(2)*Q*Q*lQ*lQ);
    H += double(1) / (Q*Q*lQ*lQ*lQ);

    return H;
  }



} // ANTL

