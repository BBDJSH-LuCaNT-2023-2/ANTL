/**
 * @file L_function_util.hpp
 * @author Leonard Nooy
 * @version $Header$
 *
 * Purpose: This file contains all of the function calls and other stuff that
 *          doesn't really belong in any other file. I know this is not a good
 *          way to do things, but hey... I'm a busy man.
 */

#include <NTL/GF2EX.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pEXFactoring.h>

#include <ANTL/common.hpp>
#include <ANTL/Arithmetic/CC.hpp>

namespace ANTL
{
  /*
   * Function: Bach_Table
   * Purpose: This function returns the A,B from the table
   *          given by Eric Bach in his paper.
   *
   * Inputs: long Q - The number of terms in the truncated Euler product
   *         double & A - A pointer to hold the value of A.
   *         double & B - A pointer to hold the value of B.
   * Outputs: NONE.
   */
  void Bach_Table (long Q, double &A, double &B);

  /*
   * Function: isInteger
   * Description: This function looks at each of the components of a complex number
   *              and checks to see if they are less that diff distance away from
   *              an integer.
   * Inputs: CC & X - The complex number to test
   *         double diff - The magnitude that the components can differ from an integer
   * Outputs: bool - 1 = X is withing diff
   *                 0 = X is not within diff, but does not say which component violated
   *                       the conditional
   */
  bool isInteger (CC < RR > &X, double diff);

  void testInteger (CC < RR > &X, RR & diff);

  /*
   * Function: Cornacchia
   * Description: This function factors a prime p=5 mod 8, such that
   *              p = a^2 + b^2, with a=-1 mod 4, b=2 mod 4, ab=2 mod 8
   * Inputs: ZZ & p - A reference to the prime to split
   *         ZZ & x - A feference to the first factor
   *         ZZ & y - A reference to the second factor
   * Outputs: void
   */
  void Cornacchia (ZZ & p, ZZ & x, ZZ & y);

  /*
   * Function: Cornacchia
   * Description: This function factors a prime p=5 mod 8, such that
   *              p = a^2 + b^2, with a=-1 mod 4, b=2 mod 4, ab=2 mod 8
   * Inputs: long & p - A reference to the prime to split
   *         long & x - A feference to the first factor
   *         long & y - A reference to the second factor
   * Outputs: void
   */
  void Cornacchia (long &p, long &x, long &y);

  /*
   * Function: Cornacchia
   * Description: This function factors a prime p=5 mod 8, such that
   *              p = a^2 + b^2, with a=-1 mod 4, b=2 mod 4, ab=2 mod 8
   * Inputs: long & p - A reference to the prime to split
   *         long & x - A feference to the first factor
   *         long & y - A reference to the second factor
   * Outputs: void
   */
  void Cornacchia (long long &p, long long &x, long long &y);

  //
  // Misc. Methods
  //

  int compute_r2(const ZZ_pX & D);
  int compute_r2(const zz_pX & D);
  int compute_r2(const ZZ_pEX & D);
  int compute_r2(const zz_pEX & D);
  int compute_r2(const GF2EX & D);
  double eval_c(const RR & D);
  double eval_U(double Q);
  double eval_G(double Q, double lQ, double U, double l);
  double eval_H(double Q, double lQ, double U, double l);

} // ANTL
