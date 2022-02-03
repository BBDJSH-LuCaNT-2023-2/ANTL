/**
 * @file L_function_FF_impl.hpp
 * @author Leonard Nooy
 * @remarks This file is to be included from L_function.h only.
 * @version $Header$
 *
 * This file contains the code that the compiler will instantiate when the
 * user attempts to create a class that is one of the NTL types: ZZ_px,
 * zz_px, ZZ_pEX, zz_PEX.
 */

namespace ANTL
{

  /***************************************************************************************************

L FUNCTIONS OVER zz_pEX, zz_pE, ZZ_pEX, ZZ_pE

  ***************************************************************************************************/

  /*
   * Function: L_function<T>::calculate_optimal_terms_ff
   * Description: This function calculates the optimal number of terms
   *              (represented as max degree of prime ideals to use for approx)
   *              to approximate L(1,X) to accuracy acc
   * Inputs: double acc - desired accuracy
   * Outputs: int - the max degree to use for the approximation
   */
  template < class T > long L_function < T >::calculate_optimal_terms_ff (double acc)
  {
    int terms = 1;
    RR err;

    err = calculate_L1_error_ff(Delta,(long) terms);
    if (info > 2)
      cout << "n = " << terms << ", acc = " << acc << ", err = " << log(err) << endl;

    while (log(err) >= acc) {
      ++terms;
      err = calculate_L1_error_ff(Delta,(long) terms);
      if (info > 2)
	cout << "n = " << terms << ", acc = " << acc << ", err = " << log(err) << endl;
    }

    return terms;
  }


  /*
   * Function: L_function<T>::calculate_L1_error_ff
   * Description: This function calculates the error between the approximate value for L1 and the
   *              real value of L1
   * Inputs: long terms - The number of terms used in the calculation of the approximation of L1
   * Outputs: RR - The error between the approximation and the real L1
   */
  template < class T > RR L_function < T >::calculate_L1_error_ff (const T & D, long Q)
  {
    ZZ q = CARDINALITY < T > ();
    ZZ g = to_ZZ (deg (Delta) >> 1);

    /* temporary variables to aid in the calculation */
    ZZ temp1num,temp1den,temp2num,temp2den;
    RR temp1RR,temp2RR;
    RR err;

    temp1num = (g << 1) + (Q % 2);
    temp2num = (g << 1) + 2;

    power(temp1den,q,Q >> 1);
    temp1den *= (Q+1);

    power(temp2den,q,(Q-1) >> 1);
    temp2den *= (Q+2);

    temp1RR = SqrRoot(to_RR(q)) - 1;
    temp2RR = temp1RR*temp1RR*temp1RR;

    err = to_RR(temp1num) / (temp1RR*to_RR(temp1den));
    err += to_RR(temp2num) / (temp2RR*to_RR(temp2den));

    return exp(err);
  }

  /*
   * Function: L_function<T>::euler_term
   * Description: Calculate the partial euler term in the Stein, Teske truncated
   *              Euler product calculation of L(s)
   * Inputs: long s - the degree of the L function
   *         long n - the degree of the polynomials to calculate over.
   * outputs: RR - the value of the euler product.
   */
  template < class T > RR L_function < T >::euler_term (long s, long n)
  {
    RR F, Fs, Ft, tempR;
    long jac;
    ZZ q, qn, LL, H;
    T P, cterm;
    long sn = 0, tn = 0;

    q = CARDINALITY < T > ();
    power (qn, q, n);
    H = (qn << 1);
    for (LL = qn; LL < H; ++LL)
      {
	get_poly_modq (P, LL, q);

	if (DetIrredTest (P))
	  {
	    jac = Chi.quadratic (P);

	    if (jac == 1)
	      ++sn;
	    if (jac == -1)
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



  /*
   * Function: L_function<T>::approximateL1
   * Description: This function calculates L(1) over a finite function field.
   *              using b truncated Euler terms.
   * Inputs: long n - the number of terms in the main truncated euler product to calculate
   * Outputs: RR - the value of L(1)
   */
  template < class T > RR L_function < T >::approximateL1 (long n)
  {
    register long i;
    ZZ q2g;			// q to the power of g + r2
    long g = deg(Delta);

    if (deg(Delta) & 1)
      g = (g-1) >> 1;
    else
      g = (g-2) >> 1;

    if (n <= nterms_used[1])
      return L1_result;

    ZZ q = CARDINALITY < T > ();

    // if this is a fresh calculation, initalize the result to
    // q^(g+r-1+r2) / [ (q-1)^(r-1) (q+1)^r2 ]

    if (nterms_used[1] == 0)
      {
	power (q2g, q, r2 + g);

	L1_result = to_RR(q2g);

	if (r2 > 0) {
	  power(q2g, q+1, r2);
	  L1_result /= to_RR(q2g);
	}

	if (!(deg(Delta) & 1)) {
	  if (test_Dcoeff(Delta,hx)) {
	    // real function field - r=2
	    L1_result *= to_RR(q) / (to_RR(q-1));
	  }
	  else {
	    // imaginary function field with deg(D) even
	    L1_result *= to_RR(2);
	  }
	}

	if (info > 2)
	  cout << "L1 init = " << L1_result << endl;
      }

    for (i = nterms_used[1] + 1; i <= n; ++i)
      {
	L1_result *= euler_term (1, i);
	if (info > 2)
	  cout << "n=" << i << ", L1 = " << L1_result << endl;
      }

    nterms_used[1] = n;

    return L1_result;
  }


  /*
   * Function: L_function<T>::approximateL1
   * Description: This function calculates L(1) over a finite function field
   *              iterativly to the given accuracy
   * Inputs: double acc - The required accuracy
   * Outputs: RR - The result of the calculation
   */
  template < class T > RR L_function < T >::approximateL1 (double acc)
  {
    long n = calculate_optimal_terms_ff(acc);
    if (info > 1)
      cout << "L1 approx with deg <= " << n << endl;
    return approximateL1(n);
  }

  /*
   * Function: L_function<T>::approximateL1
   * Description: This function calculates L(s) over a finite function field
   *              using the given number of terms.
   * Inputs: long s - The value of s in the L function
   *         long terms - The number of terms to use in the summation
   * Outputs: RR - The result of the calculation
   */
  template < class T > RR L_function < T >::approximateL (long s, long terms)
  {
    register long i;
    ZZ q2g;			// q to the power of g/2
    RR result;

    ZZ q = CARDINALITY < T > ();
    ZZ qs;			// q to the power of s

    power (qs, q, s);
    power (q2g, qs, deg (Delta) / 2);

    result = to_RR (q2g) / to_RR (qs - 1);

    for (i = 1; i <= terms; ++i)
      {
	result *= euler_term (s, i);
      }

    return result;
  }

  /*
   * approximate the value of L(s)
   *
   */
  template < class T > RR L_function < T >::approximate (long s, long n)
  {
    switch (s)
      {
      case 1:
	return approximateL1 (n);

      default:
	return approximateL (s, n);
      }
  }

  template < class T > RR L_function < T >::approximate (long s, double acc)
  {
    switch (s)
      {
      case 1:
	return approximateL1 (acc);

      default:
	return approximateL ((long) 6);	// the number of terms to use here should really be calculated somewhere.
      }
  }

} // ANTL
