/**
 * @file L_function_GF2EX.cpp
 * @author Leonard Nooy
 * @version $Header$
 *
 * This file contains the GF2EX specialization for the L_function class
 */

#include <ANTL/L_function/L_function.hpp>

namespace ANTL
{

  /*
  //
  // L_function<GF2EX>::euler_term1()
  //
  // Task:
  //      returns the euler term of L(1) for degree 1
  //

  template <> RR quadratic_order < GF2EX >::euler_term1 (long ss) const
  {
  RR F, Fs, Ft, tempR;
  long jac;
  ZZ q, LL;
  long s = 0, t = 0;

  GF2EX P, x, cterm;

  SetX (x);

  q = CARDINALITY<GF2EX>();
  for (LL = 0; LL < q; ++LL)
  {
  get_poly_modq (cterm, LL, q);
  P = x + cterm;

  jac = Jacobi (hx, Delta, P);

  if (jac == 1)
  ++s;
  if (jac == -1)
  ++t;
  }

  // evaluate Euler term
  ZZ qn;

  power (qn, qn, ss);

  tempR = to_RR (q) / to_RR (q - 1);
  power (Fs, tempR, s);

  tempR = to_RR (q) / to_RR (q + 1);
  power (Ft, tempR, t);

  F = Fs * Ft;

  return F;
  }



  //
  // L_function<GF2EX>::euler_term2()
  //
  // Task:
  //      returns the euler term of L(1) for degree 2
  //

  template <> RR quadratic_order < GF2EX >::euler_term2 (long ss) const
  {
  RR F, Fs, Ft, tempR;
  long jac;
  ZZ q, qn, LL, H;
  GF2EX P, x, cterm;
  long s = 0, t = 0;

  SetX (x);

  q = CARDINALITY<GF2EX>();
  qn = q * q;
  H = (qn << 1);
  for (LL = qn; LL < H; ++LL)
  {
  get_poly_modq (P, LL, q);

  if (DetIrredTest (P))
  {
  jac = Jacobi (hx, Delta, P);

  if (jac == 1)
  ++s;
  if (jac == -1)
  ++t;
  }
  }

  // evaluate Euler term
  power (qn, qn, ss);

  tempR = to_RR (qn) / to_RR (qn - 1);
  power (Fs, tempR, s);

  tempR = to_RR (qn) / to_RR (qn + 1);
  power (Ft, tempR, t);

  F = Fs * Ft;

  return F;
  }
  */


  /*
   * Function: euler_term
   * Purpose: Calculate the Euler term of degree n, for the calculation
   *          of L(1).
   * Inputs: long s - For L(1) this will be set to 1,
   *                  For L(2) this will be set to 2, etc...
   *         long n - the degree of the monic prime polynomials to evaluate
   *                  this Euler term over.
   */

  template <> RR L_function < GF2EX >::euler_term (long s, long n)
  {
    RR F, Fs, Ft, tempR;
    long art;
    ZZ q, qn, LL, H;
    GF2EX P, cterm;
    long sn = 0, tn = 0;

    /*
      if (n == 1)
      return euler_term1 (ss);

      if (n == 2)
      return euler_term2 (ss);
    */

    q = CARDINALITY<GF2EX>();
    power (qn, q, n);

    H = (qn << 1);
    for (LL = qn; LL < H; ++LL)
      {
	get_poly_modq (P, LL, q);

	if (DetIrredTest (P))
	  {
	    art = Chi.quadratic (P);

	    if (art == 1)
	      ++sn;
	    if (art == -1)
	      ++tn;
	  }
      }

    // evaluate Euler term
    power (qn, qn, s);

    tempR = to_RR (qn) / to_RR (qn - 1);
    power (Fs, tempR, sn);

    tempR = to_RR (qn) / to_RR (qn + 1);
    power (Ft, tempR, tn);

    F = Fs * Ft;

    return F;
  }

} // ANTL

