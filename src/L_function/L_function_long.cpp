/**
 * @file L_function_long.cpp
 * @author Leonard Nooy
 * @version $Header$
 *
 * This file contains the long specialization for the L_function class.
 */

#include <ANTL/L_function/L_function.hpp>

#ifdef DEBUG
#include <typeinfo>
#endif

using namespace std;

namespace ANTL
{

  /*******************************************************************************************

        INITALIZATION FUNCTIONS

  *******************************************************************************************/

  /*
   * Function: L_function<long>::init
   * Description: This function must be called before using the rest of the L_function class
   * Inputs: const long & x - The conductor of the calculation to preform
   *         int m - The mode in which the class should run (This only effects the L0 calculations )
   * Outputs: none.
   */

  template <> void L_function < long >::init (const long &x, int m)
  {
    if (m == QUADRATIC_MODE || m == QUARTIC_MODE)
      mode = m;
    else
      {
        cout << "Init: Unrecognized Mode, Defaulting to Quadratic Mode." << endl;
        mode = QUADRATIC_MODE;
      }

    Delta = x;
    memset (nterms_used, 0, sizeof (long) * 3);
    Chi.init(Delta,mode);

    /* if the mode if quartic set up the a,b for the L0 calculation */
    if (mode == QUARTIC_MODE)
      {
        Cornacchia (Delta, a, b);
      }
  }


  /*******************************************************************************************

        L0 FUNCTIONS

  *******************************************************************************************/

  /************************************* QUARTIC *********************************************/

  /*
   * Function: template <class TYPE> approximateL0_Quartic_impl
   * Description: This function calculates L0 in the Imaginary Quartic case. It uses techniques
   *              divised by Stephean Louboutin.
   *              Generic template version by Rennie deGraaf - computes using the precision of
   *              floating-point type TYPE.  This function is a friend of the L_function class.
   * Inputs:      L_function<long>& lfunc - the L-function to compute
   *              long terms - The number of terms to use in the summation.
   * Outputs: CC - The result of the calculation
   */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<long>& lfunc, long terms)
  {
    long i;		// loop var
    CC<TYPE>T1, temp3;	// used in the calulation of the two multipliers
    TYPE T2, temp1, temp2;
    TYPE f, h, en;		// used to determine en inductively
    CC<TYPE> S1, S2, S3;	// these represent the three sums in the calculation
    CC<TYPE> CHI, CHI_HAT;	// the values from the chi character

    /* Calculate the multipliers for the two sums:
     * T1 = sqrt(p + a*sqrt(p) / 2) + i * b / |b| * sqrt(p - a*p / 2)
     * T2 = 1/sqrt(p)
     */
#ifdef DEBUG
    cout.precision(10);
    cout << "Initial state: Delta = " << lfunc.Delta << " PI = " << lfunc.PI
	 << " a = " << lfunc.a << " b = " << lfunc.b << endl;
    cout << "Instance: " << typeid(TYPE).name() << endl;
#endif
    T2 = sqrt(to<TYPE>(lfunc.Delta));
#ifdef DEBUG
    cout << "Init 1: " << to<TYPE>(lfunc.Delta) << ", " << sqrt(to<TYPE>(lfunc.Delta)) << endl;
    cout << "Delta: " << typeid(lfunc.Delta).name() << " to<TYPE>(Delta): " << typeid(to<TYPE>(lfunc.Delta)).name() << endl;
    cout << "(RR)Delta: " << to<RR>(lfunc.Delta) << " sqrt((RR)Delta): " << sqrt(to<RR>(lfunc.Delta)) << endl;
    cout << "(float)Delta: " << to<float>(lfunc.Delta) << " sqrt((float)Delta): " << sqrt(to<float>(lfunc.Delta)) << endl;
    cout << "(double)Delta: " << to<double>(lfunc.Delta) << " sqrt((double)Delta): " << sqrt(to<double>(lfunc.Delta)) << endl;
#endif
    temp1 = sqrt((to<TYPE>(lfunc.Delta) + to<TYPE>(lfunc.a) * T2) / to<TYPE>(2.0));
    temp2 = sqrt((to<TYPE>(lfunc.Delta) - to<TYPE>(lfunc.a) * T2) / to<TYPE>(2.0)) *
      (to<TYPE>(lfunc.b) / abs(to<TYPE>(lfunc.b)));
    T1.assign(temp1, temp2);
#ifdef DEBUG
    cout << "Init 4: T2 = " << T2 << " temp1 = " << temp1 << " temp2 = " << temp2
	 << " T1 = " << T1 << endl;
#endif
    temp3.assign(to<TYPE>(0), to<TYPE>(lfunc.PI));
    divide(T1, T1, temp3);

    T2 = to<TYPE>(1) / T2;

    /* setup the variables for the calculation of en */
    f = exp(to<TYPE>(-lfunc.PI) / to<TYPE>(lfunc.Delta));
    h = f * f;
    en = f;

    S1.clear();
    S2.clear();
    S3.clear();

#ifdef EXTRA
    cout << "a % 3 = " << (lfunc.a % 3) << endl;
    cout << "b % 5 = " << (lfunc.b % 5) << endl;
#endif

    long pcount = 1;

    for (i = 1; i < terms; i++)
      {
#ifdef DEBUG
	cout << "Cycle " << i << ": T1 = " << T1 << " T2 = " << T2 << " S1 = " << S1
	     << " S2 = " << S2 << " S3 = " << S3 << " en = " << en << " f = " << f
	     << " h = " << h << endl;
#endif
        // use en in the first sum,
        if (pcount == lfunc.Delta)
	  {
            pcount = 0;
            CHI.clear();
            CHI_HAT.clear();
	  }
        else
	  {
#ifdef USE_CHI2
            CHI = to< CC<TYPE> >(lfunc.Chi.quartic(i));
#else
            CHI = to< CC<TYPE> >(lfunc.Chi.quartic_table(i));
#endif
            //CHI_HAT.conjugate(CHI);
            conjugate(CHI_HAT, CHI);
	  }
        ++pcount;

#ifdef EXTRA
        if (i < 100 && ProbPrime (i))
	  cout << "i=" << i << ", X = " << CHI << endl;
#endif

        temp1 = en / to<TYPE>(i);

        multiply(CHI_HAT, CHI_HAT, temp1);
        add(S1, S1, CHI_HAT);
        add(S3, S3, CHI);		// S3

        temp1 = en;

        // calculate en+1 for the S2 sum
        f = h * f;
        en = en * f;

        temp1 += en;
        multiply(CHI, S3, temp1);
        add(S2, S2, CHI);
      }

    /* we now have the values of S1, S2.
     * so Calculate T1*S1 + T2*S2, and T1*S1 - T2*S2, and see which one is an integer.
     */

    // L1 = T2 S2 + T1 S1, L2 = T2 S2 - T1 S1
    CC <TYPE> L_1, L_2;

    multiply(L_1, S2, T2);
    L_2 = L_1;
    multiply(temp3, T1, S1);
    add(L_1, L_1, temp3);
    subtract(L_2, L_2, temp3);

    CC<RR> L_1_RR = to< CC<RR> >(L_1);
    CC<RR> L_2_RR = to< CC<RR> >(L_2);

#ifdef EXTRA
    cout << "\np=" << lfunc.Delta << " (n=" << terms << "), L1 = " << L_1 << ", L2 = "
	 << L_2 << endl;
#endif

    RR diff1, diff2;

    testInteger(L_1_RR, diff1);
    testInteger(L_2_RR, diff2);

    if (diff1 > 0.015 && diff2 > 0.015)
      {
        cerr << "Error: L1 and L2 both are not close to being integers! Please run with more terms and/or more precision."
             << endl;
        cout << "p = " << lfunc.Delta << endl;
        lfunc.epsilon = 0;
      }

    if (diff1 < diff2)
      {
        lfunc.L0_result = L_1_RR;
        lfunc.nterms_used[L0_REF] = terms;
        lfunc.epsilon = 1;
      }
    else
      {
        lfunc.L0_result = L_2_RR;
        lfunc.nterms_used[L0_REF] = terms;
        lfunc.epsilon = -1;
      }

    lfunc.L0_result.assign(round(lfunc.L0_result.real()), round(lfunc.L0_result.imaginary()));

#ifdef EXTRA
    cout << "L0 = " << lfunc.L0_result << endl;
#endif

    return lfunc.L0_result;
  }


  /*
   * Function: L_function<long>::approximateL0_Quartic
   * Description: This function is simply a wrapper for the calculate_optimal_n and approximateL0(long n)
   * Inputs: double acc - The required accuracy
   * Outputs: CC - The result of the calculation
   */
  template <> CC < RR > L_function < long >::approximateL0_Quartic (double acc)
  {
    return approximateL0_Quartic_default(acc);
  }

  /*********************************** IMAGINARY NUMBER FIELDS ******************************/

  /*
   * Function: template <class TYPE> approximateL0_ImaginaryNumberField_impl
   * Description: This function calculates L0 in the Imaginary Quadratic case. It uses techniques
   *              divised by Stephean Louboutin.
   *              Generic template version by Rennie deGraaf - computes using the precision of
   *              floating-point type TYPE.  This function is a friend of the L_function class.
   * Inputs:      L_function<long>& lfunc - the L-function to compute
   *              long terms - The number of terms to use in the summation.
   * Outputs: RR - The result of the calculation
   */

  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<long>& lfunc, long terms)
  {
    long int i;			// loop var
    TYPE T1, T2;			// used in the calulation of the two multipliers
    TYPE temp1, temp2;
    TYPE f, h, en;		// used to determine en inductively
    TYPE S1, S2, S3;		// these represent the three sums in the calculation
    long CHI;			// the values from the chi character
    long dk;

    /* Calculate the multipliers for the two sums:
     * T1 = sqrt(Delta) / PI
     * T2 = sqrt(PI/Delta)
     */
    dk = labs(lfunc.Delta);
    T1 = sqrt(to<TYPE>(dk)) / to<TYPE>(lfunc.PI);
    T2 = to<TYPE>(1) / sqrt(to<TYPE>(dk));

    /* setup the variables for the calculation of en */
    f = exp(to<TYPE>(-lfunc.PI) / to<TYPE>(dk));
    h = f * f;
    en = f;

    clear(S1);
    clear(S2);
    clear(S3);
    clear(temp1);
    clear(temp2);

    for (i = 1; i < terms; i++)
      {
        CHI = lfunc.Chi.quadratic(i);

        // use en in the first sum,
        temp1 = en / to<TYPE>(i);
        temp1 *= CHI;
        S1 += temp1;

        S3 += CHI;

        temp1 = en;		// S2
        // calculate en+1 for the S2 sum
        f = h * f;
        en = en * f;

        temp1 += en;
        temp2 = S3 * temp1;
        S2 += temp2;
      }

    lfunc.L0_result.assign (to_RR ((S1 * T1) + (S2 * T2)), to_RR (0));
    lfunc.nterms_used[L0_REF] = terms;

    return lfunc.L0_result.real ();
  }



  /*
   * Function: L_function<long>::approximateL0_ImaginaryNumberField
   * Description: This function is simply a wrapper for the calculate_optimal_n and approximateL0(long n)
   * Inputs: double acc - The required accuracy
   * Outputs: CC - The result of the calculation
   */
  template <> RR L_function < long >::approximateL0_ImaginaryNumberField (double acc)
  {
    return approximateL0_ImaginaryNumberField_default(acc);
  }

  /*************************************** REAL NUMBER FIELDS *****************************************/

  /*
   * Function: template <class TYPE> approximateL0_RealNumberField_impl
   * Description: This function calculates L0 in the Real Quartic case. It uses techniques
   *              divised by Stephean Louboutin.
   *              Generic template version by Rennie deGraaf - computes using the precision of
   *              floating-point type TYPE.  This function is a friend of the L_function class.
   * Inputs:      L_function<long>& lfunc - the L-function to compute
   *              long terms - The number of terms to use in the summation.
   * Outputs: RR - The result of the calculation
   */
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<long>& lfunc, long terms)
  {
    long int i;			// loop var
    TYPE temp1;			// used in the calulation of the two multipliers
    TYPE f, h, en, en_1;		// used to determine en inductively
    TYPE S1, S2, S3, S4;		// these represent the three sums in the calculation
    long CHI;			// the values from the chi character

    /* setup the variables for the calculation of en */
    f = exp(to<TYPE>(-lfunc.PI) / to<TYPE>(lfunc.Delta));
    h = f * f;
    en = f;

    clear(S1);
    clear(S2);
    clear(S3);
    clear(S4);
    clear(temp1);

    for (i = 1; i < terms; i++)
      {
        CHI = lfunc.Chi.quadratic(i);

        // calculate en, and en + 1
        en_1 = en;
        f = h * f;
        en = en * f;

        // calculate S3, S4
        S4 += to<TYPE>(CHI) / to<TYPE>(i);
        S3 += (en + en_1) * S4;

        // calculate S2, S1
        temp1 = ((en_1 - to<TYPE>(1)) / to<TYPE>(i)) + ((en - to<TYPE>(1)) / to<TYPE>(i + 1));
        temp1 += to<TYPE>(2) * log (to<TYPE>(i + 1) / to<TYPE>(i));
        S2 += CHI;
        S1 += temp1 * S2;
      }

    lfunc.L0_result.assign(to_RR(S1 + S3) / to_RR(2), to_RR(0));
    lfunc.nterms_used[L0_REF] = terms;

    return lfunc.L0_result.real ();
  }



  /*
   * Function: L_function<long>::approximateL0_RealNumberField
   * Description: This function is simply a wrapper for the calculate_optimal_n and approximateL0(long n)
   * Inputs: double acc - The required accuracy
   * Outputs: CC - The result of the calculation
   */
  template <> RR L_function < long >::approximateL0_RealNumberField (double acc)
  {
    return approximateL0_RealNumberField_default(acc);
  }


  /* for the time being this is function only calculates using the quad_float base type */
  template <> CC < RR > L_function < long >::approximateL0 (long terms)
  {
    return approximateL0_default(terms);
  }


  template <> CC < RR > L_function < long >::approximateL0 (double acc)
  {
    return approximateL0_default(acc);
  }

  /*******************************************************************************************

        L1 FUNCTIONS

  *******************************************************************************************/

  /*
   * Function: template <class TYPE> approximateL1_Quartic_impl
   * Description: This function calculates the value of L1 using eric bach's method
   *              of truncated euler terms. It used doubles for internal calculations
   *              Generic template version by Rennie deGraaf - computes using the precision of
   *              floating-point type TYPE.  This function is a friend of the L_function class.
   * Inputs:      L_function<long>& lfunc - the L-function to compute
   *              long terms - The number of terms to use in the summation.
   * Outputs: RR - The result of the calculation
   */
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<long>& lfunc, double acc)
  {
    long P, kron;
    TYPE E, nE, P2, term;
    RR Epref;
    CC <float> X;

    Epref = to_RR(lfunc.Delta) / (to_RR(2) * lfunc.PI * lfunc.PI);
    E = to<TYPE>(0.8);
#ifdef EXTRA
    cout << "P = 2 - " << Epref*to_RR(E) << endl;
#endif

    // compute partial product  p < Q
    lfunc.primes.reset(3);
    P = lfunc.primes.next ();		// primes start at 3
    while (1)
      {
        P2 = to<TYPE>(P) * to<TYPE>(P);

        kron = Kronecker(lfunc.Delta,P);
        if (kron == -1) {
	  term = P2 / (P2 + 1);
#ifdef EXTRA
	  if (P < 100) {
	    cout << "\ni = " << P << ", X = i,-i" << endl;
	    cout << "L = " << E / (to_RR(lfunc.Delta) / (to_RR(2)*lfunc.PI*lfunc.PI)) << endl;
	  }
#endif
        }
        else if (kron == 1)
	  {
            X = lfunc.Chi.quartic(P);
            if (X.is_one())
	      term = P2 / (P2 - (P << 1) + 1);
            else
	      term = P2 / (P2 + (P << 1) + 1);
#ifdef EXTRA
            if (P < 100) {
	      cout << "\ni = " << P << ", X = " << X << endl;
	      cout << "L = " << E / (to_RR(lfunc.Delta) / (to_RR(2)*lfunc.PI*lfunc.PI)) << endl;
            }
            if (X.imaginary() != 0)
	      cerr << "ERROR!  complex character!" << endl;
#endif
	  }
        else
	  term = to<TYPE>(1);

        nE = E*term;
        if (kron == 1 && abs(Epref*to_RR(nE-E)) < acc && abs(Epref*to_RR(E) - round(Epref*to_RR(E))) < 0.01)
	  break;

        E = nE;
#ifdef EXTRA
        cout << "P = " << P << " - " << E << ", L = " << to_double(sqrt(E / (to_RR(lfunc.Delta) / (to_RR(2)*lfunc.PI*lfunc.PI)))) <<endl;
#endif
        P = lfunc.primes.next ();
      }

    lfunc.nterms_used[L1_REF] = P;

    return round(Epref*to_RR(E));
  }


  /*
   * Function: template <class TYPE> approximateL1Bach_Quartic_impl
   * Description: This function calculates the value of L1 using eric bach's method
   *              of truncated euler terms. It used doubles for internal calculations
   *              Generic template version by Rennie deGraaf - computes using the precision of
   *              floating-point type TYPE.  This function is a friend of the L_function class.
   * Inputs:      L_function<long>& lfunc - the L-function to compute
   *              long terms - The number of terms to use in the summation.
   * Outputs: RR - The result of the calculation
   */
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<long>& lfunc, long terms)
  {
    long Q, Q2, P, kron;
    TYPE C, E, wt, term, dP;
    long i;
    CC<float> X;

    Q = terms;

    // compute weight
    Q2 = Q << 1;			// Q2 = Q*2;
    clear (C);
    for (i = Q; i <= Q2 - 1; ++i)
      C += to<TYPE>(i) * log(to<TYPE>(i));

    // compute partial product  p < Q
    E = log(to<TYPE>(0.8));     // assuming p = 5 mod 8
    lfunc.primes.reset(3);
    P = lfunc.primes.next();		// primes start at 2
    while (P < Q)
      {
        dP = to<TYPE>(P);
        kron = Kronecker(lfunc.Delta,P);
        if (kron == -1) {
	  term = dP / sqrt(dP*dP+1);
#ifdef EXTRA
	  if (P < 500) {
	    cout << "\ni = " << P << ", X = i,-i" << endl;
	    cout << "L = " << exp(to_RR(E)) << endl;
	  }
#endif
        }
        else if (kron == 1)
	  {
            X = lfunc.Chi.quartic(P);
            if (X.is_one())
	      term = dP / (dP - 1);
            else
	      term = dP / (dP + 1);
#ifdef EXTRA
            if (P < 500)
	      cout << "i = " << P << ", X = " << X << endl;
#endif
	  }
        else
	  set(term);

        E += log(term);
        P = lfunc.primes.next();
      }

    // computed weighted partial products for Q < p < 2Q
    set(wt);
    for (i = Q; i <= P; ++i)
      wt -= to<TYPE>(i) * log(to<TYPE>(i)) / C;

    while (P < Q2)
      {
        dP = to<TYPE>(P);
        kron = Kronecker(lfunc.Delta,P);
        if (kron == -1) {
	  term = dP / sqrt(dP*dP+1);
#ifdef EXTRA
	  if (P < 500)
            cout << "\ni = " << P << ", X = i,-i" << endl;
#endif
        }
        else if (kron == 1)
	  {
            X = lfunc.Chi.quartic(P);
            if (X.is_one())
	      term = dP / (dP - 1);
            else
	      term = dP / (dP + 1);
#ifdef EXTRA
            if (P < 500)
	      cout << "i = " << P << ", X = " << X << endl;
#endif
	  }
        else
	  set(term);

        E += wt * log(term);
        P = lfunc.primes.next();
        wt -= (to<TYPE>(P - 1) * log(to<TYPE>(P - 1)) + to<TYPE>(P) * log(to<TYPE>(P))) / C;
      }

    lfunc.nterms_used[L1_REF] = Q;

    return (exp(to_RR(E)));
  }



/*
* Function: template <class TYPE> approximateL1_impl
* Description: This function calculates the value of L1 using eric bach's method
*              of truncated euler terms. It used doubles for internal calculations
*              Generic template version by Rennie deGraaf - computes using the precision of
*              floating-point type TYPE.  This function is a friend of the L_function class.
* Inputs:      L_function<long>& lfunc - the L-function to compute
*              long terms - The number of terms to use in the summation.
* Outputs: RR - The result of the calculation
*/
template <class TYPE> RR approximateL1_impl(L_function<long> &lfunc, long terms) {
  long Q, Q2, P, kron;
  TYPE C, E, wt;
  long i;

  Q = terms;

  // compute weight
  Q2 = Q << 1; // Q2 = Q*2;
  clear(C);
  for (i = Q; i <= Q2 - 1; ++i)
    C += to<TYPE>(i) * log(to<TYPE>(i));

  // compute partial product  p < Q
  clear(E);
  lfunc.primes.reset(2);
  P = lfunc.primes.next(); // primes start at 2
  while (P < Q) {
    kron = lfunc.Chi.quadratic(P);
    E += log(to<TYPE>(P) / to<TYPE>(P - kron));
    P = lfunc.primes.next();
  }

  // computed weighted partial products for Q < p < 2Q
  set(wt);
  for (i = Q; i <= P; ++i)
    wt -= to<TYPE>(i) * log(to<TYPE>(i)) / C;

  while (P < Q2) {
    kron = lfunc.Chi.quadratic(P);
    E += wt * log(to<TYPE>(P) / to<TYPE>(P - kron));
    P = lfunc.primes.next();
    wt -= (to<TYPE>(P - 1) * log(to<TYPE>(P - 1)) +
           to<TYPE>(P) * log(to<TYPE>(P))) /
          C;
  }

  lfunc.nterms_used[L1_REF] = Q;

  return (exp(to_RR(E)));
}



  /*
   * Function: L_function<long>::approximateL1
   * Description: This function calculates the L(1,X) function over the finite
   *              number field of quadratic order. The descriminant is
   *              represented by a long, and thus has no restriction.
   *              This function used Eric Bach's method of weighted truncated
   *              Euler products.
   * Inputs: long n - the number of term in the truncated Euler product
   * Output: RR - the result of the product.
   */
  template <> RR L_function < long >::approximateL1 (long n)
  {
    /* if we have precalculated the L1 function to this many terms return that value */
    if (n <= nterms_used[L1_REF])
      {
        return L1_result;
      }

    long Q;
    if (global_flag)
      Q = nterms_used[L1_REF] = global_num_terms;
    else
      Q = nterms_used[L1_REF] = n;

    if (mode == QUADRATIC_MODE) {
      long prec = calculate_precision (Delta);
      if (prec <= 53) {
        L1_result = approximateL1_impl<double>(*this, Q);
      }
      else if (prec <= 104) {
        L1_result = approximateL1_impl<quad_float>(*this, Q);
      }
      else {
        RR::SetPrecision (prec);
        L1_result = approximateL1_impl<RR>(*this, Q);
      }
    }
    else if (mode == QUARTIC_MODE) {
        long prec = calculate_precision(Delta);
        RR::SetPrecision(NumBits(Delta) << 1);
        RR::SetOutputPrecision(NumBits(Delta) << 1);
        ComputePi(PI);

        if (prec <= 53)
	  {
            L1_result = approximateL1Bach_Quartic_impl<double>(*this, n);
	  }
        else if (prec <= 104)
	  {
            L1_result = approximateL1Bach_Quartic_impl<quad_float>(*this, n);
	  }
        else
	  {
	    //            RR::SetPrecision (prec);
            L1_result = approximateL1Bach_Quartic_impl<RR>(*this, n);
	  }
    }

    return L1_result;
  }

  /*
   * Function: L_function<long>::approximateL1
   * Description: This function calls calculate_optimal_terms to determine how many terms
   *              are required to keep the error below the given accuracy, and also calculates
   *              the amount of roundoff error that is acceptable. Based on that it chooses
   *              an appropriate precision specialization of the approximateL1 to do the
   *              actual calculation.
   * Inputs: double acc - The required accuracy of the calculation
   * Outputs: RR - The result of the calculation
   */
  template <> RR L_function < long >::approximateL1 (double acc)
  {
    long Q;
    long prec;

    if (mode == QUADRATIC_MODE)
      {
        if (global_flag)
	  {
            Q = global_num_terms;
            prec = global_prec;

            if (prec <= 53)
	      {
                L1_result = approximateL1_impl<double>(*this, Q);
	      }
            else if (prec <= 104)
	      {
                L1_result = approximateL1_impl<quad_float>(*this, Q);
	      }
            else
	      {
                /* there is no need to set the RR precision because it was set when the global
                   precision was setup */
                L1_result = approximateL1_impl<RR>(*this, Q);
	      }
	  }
        else
	  {
            Q = calculate_optimal_terms_new (Delta, acc);
            if (Q <= nterms_used[L1_REF])
	      {			/* if the L0 function has been calculated before, return that value */
                return L1_result;
	      }

            prec = calculate_precision (Delta);
            if (prec <= 53)
	      {
                L1_result = approximateL1_impl<double>(*this, Q);
	      }
            else if (prec <= 104)
	      {
                L1_result = approximateL1_impl<quad_float>(*this, Q);
	      }
            else
	      {
                RR::SetPrecision (prec);
                L1_result = approximateL1_impl<RR>(*this, Q);
	      }
	  }
      }
    else if (mode == QUARTIC_MODE)
      {
        prec = calculate_precision (Delta);
        RR::SetPrecision(NumBits(Delta) << 1);
        RR::SetOutputPrecision(NumBits(Delta) << 1);
        ComputePi(PI);
        prec = 200;

        if (prec <= 53)
	  {
            L1_result = approximateL1_Quartic_impl<double>(*this, acc);
	  }
        else if (prec <= 104)
	  {
	    L1_result = approximateL1_Quartic_impl<quad_float>(*this, acc);
	  }
        else
	  {
	    //	          RR::SetPrecision (prec);
            L1_result = approximateL1_Quartic_impl<RR>(*this, acc);
	  }
      }

    return L1_result;
  }

  /********************************************************************************************

         L2 FUNCTIONS

  ********************************************************************************************/

  /*
   * Function: L_function<long>::approximateL2
   * Description: This function is similar to the L1 functions, except it simply does
   *              a straight truncated product to calculate L2
   * Inputs: long terms - The number of terms in the summation
   * Outpus: RR - The result of the calculation
   */
  template <> RR L_function < long >::approximateL2 (long terms)
  {
    long kron, P;
    double E;

    if (terms <= nterms_used[L2_REF])
      {
	return L2_result;
      }

    // compute the partial product of primes less than terms
    E = to<double>(1);
    primes.reset(2);
    P = primes.next ();
    while (P < terms)
      {
	kron = Chi.quadratic(P);
	P *= P;			// P squared since s = 2
	E *= std::log (to<double>(P) / to<double>(P - kron));
	P = primes.next ();
      }

    nterms_used[L2_REF] = terms;
    L2_result = to_RR (E);
    return L2_result;
  }


  /********************************************************************************************

         GENERAL L FUNCTIONS

  ********************************************************************************************/

  /*
   * Function: L_function<long>::approximateL
   * Description: This function will calculate L(s,x) to the desired number of terms using
   *              a straight truncated product
   * Inputs: long s - The value of s (must be != 0)
   * Outputs: long terms - The number of terms to carry the calculation out to.
   */
  template <> RR L_function < long >::approximateL (long s, long terms)
  {
    long i;
    long kron, P;
    double E, Ps;

    primes.reset(2);
    P = primes.next ();
    // compute partial product for LQ <= p < QQ
    set (E);
    while (P < terms)
      {
	kron = Chi.quadratic(P);
	Ps = to<double>(P);
	for (i = 1; i < s; ++i)	// raise P to the power of s
	  Ps *= P;
	E += std::log (Ps / (Ps - kron));
	P = primes.next ();
      }

    return exp (to_RR (E));
  }

  /*
   * Function: L_function<long>::approximate
   * Description: This function is simply a wrapper for the rest of the L functions
   * Inputs: long s - The value of s in the L calculation
   * Outputs: long n - The number of terms to carry the calculation out to
   */
  template <> RR L_function < long >::approximate (long s, long n)
  {
    switch (s)
      {
      case 0:
	if (mode == QUARTIC_MODE) {
	  cout <<
	    "Please call the approximateL0 function explicitly as it returns Complex Values"
	       << endl;
	  return to_RR (0.0);
	}
	else
	  return approximateL0 (n).real();
      case 1:
	return approximateL1 (n);

      case 2:
	return approximateL2 (n);

      default:
	return approximateL (s, n);
      }
  }


  /*
   * Function: L_function<long>::approximate
   * Description: This function is simply a wrapper for the rest of the L functions
   * Inputs: long s - The value of s in the L calculation
   * Outputs: double acc - The desired accuracy of the calculation
   */
  template <> RR L_function < long >::approximate (long s, double acc)
  {
    switch (s)
      {
      case 0:
	if (mode == QUARTIC_MODE) {
	  cout <<
	    "Please call the approximateL0 function explicitly as it returns Complex Values"
	       << endl;
	  return to_RR (0.0);
	}
	else
	  return approximateL0 (acc).real();
      case 1:
	return approximateL1 (acc);

      default:
	return approximateL (s, (long) 10000);	// the value passed to the function should really be calculated
	// instead of just defaulting to 10000
      }
  }



  /******************************************************************************************

        TABLE DRIVEN FUNCTIONS

  ******************************************************************************************/

  /****************************** L0 FUNCTIONS ***********************************************/

  template <>
  RR L_function < long >::approximateL0_RealNumberField_table(double acc)
  {
    RR L0;
    long n = calculate_optimal_n_RealNumberField (Delta, acc);

    if (info > 2)
      cout << "\nStarting approxL0_real_table:  n = " << n << ", table_size = " << table_size << endl;

    // make sure the tables are created and ready to go
    if (n > table_size && table_size < MAX_L0_TABLE) {
      cerr << "ERROR (approximateL0_RealNumberField_table):  table_size set incorrectly" << endl;
      exit(1);
    }

    // make sure tables are set up
    if (!chi_table)
      {
	cerr <<
	  "Please call create_L0_tables before using the table driven functions."
	     << endl;
	exit(1);
      }

    if (global_prec == 53) {
      if (info > 2)
	cout << "Using double" << endl;

      L0 = approximateL0_RealNumberField_table_impl<double>(*this,n);
    }
    else {
      if (info > 2)
	cout << "Using quad_float" << endl;

      L0 = approximateL0_RealNumberField_table_impl<quad_float>(*this,n);
    }

    nterms_used[L0_REF] = n;

    if (info > 2)
      cout << "Done approxL0_real_table - L0 = " << L0 << endl;

    return L0;
  }



  template <class TYPE> RR approximateL0_RealNumberField_table_impl(L_function<long>& lfunc, long n)
  {
    long int i,Di;     // loop var
    //  TYPE T1, T2;	      // used in the calulation of the two multipliers
    TYPE temp1;	      // used in the calulation of the two multipliers
    TYPE f, h, en, en_1;        // used to determine en inductively
    TYPE S1, S2, S3, S4;        // the three sums in the calculation
    TYPE d2,d1,CHI,CHI_n;
    long D2,D8;

    // setup the variables for the calculation of en
    f = exp (to<TYPE>(-lfunc.PI) / to<TYPE>(lfunc.Delta));
    h = f*f;
    en = f;

    clear (S1);
    clear (S2);
    clear (S3);
    clear (S4);
    clear (temp1);

    d1 = to<TYPE>(1);
    d2 = to<TYPE>(2);

    D2 = (lfunc.Delta & 1);
    D8 = (lfunc.Delta & 7);

    for (i = 1; i < n; i++)
      {
	if (!(i & 1) && !D2) // D and p are both even
	  CHI = CHI_n = to<TYPE>(0);
	else
	  {
	    bool neg = false;
	    long t = i;

	    while (!(t & 1))
	      {
		t >>= 1;
		neg = !neg;
	      }

	    if (t > 1 && t < lfunc.table_size) {
	      Di = lfunc.Delta % t;
	      CHI = to<TYPE>(lfunc.chi_table[(t-1) >> 1][Di]);
	    }
	    else
	      CHI = lfunc.Chi.quadratic_long(t);

	    if (neg && (D8 == 5))
	      CHI = -CHI;

	    CHI_n = CHI/to<TYPE>(i);
	  }


	// calculate en, and en + 1
	en_1 = en;
	f *= h;
	en *= f;

	// calculate S3, S4
	S4 += CHI_n;
	S3 += (en + en_1) * S4;

	// calculate S2, S1
	temp1 = ((en_1 - d1) / to<TYPE>(i)) + ((en - d1) / to<TYPE>(i + 1));
	temp1 += d2 * log (to<TYPE>(i + 1) / to<TYPE>(i));
	S2 += CHI;
	S1 += temp1 * S2;
      }

    lfunc.L0_result.assign (to<RR>((S1 + S3) / d2), to<RR> (0));

    return lfunc.L0_result.real ();
  }



  template <>
  RR L_function < long >::approximateL0_ImaginaryNumberField_table(double acc)
  {
    RR L0;
    long n = calculate_optimal_n_ImaginaryNumberField (Delta, acc);

    if (info > 2)
      cout << "\nStarting approxL0_imag_table:  n = " << n << ", table_size = " << table_size << endl;

    // make sure the tables are created and ready to go
    if (n > table_size && table_size < MAX_L0_TABLE) {
      cerr << "ERROR (approximateL0_ImaginaryNumberField_table):  table_size set incorrectly" << endl;
      exit(1);
    }

    // make sure tables are set up
    if (!chi_table)
      {
	cerr <<
	  "Please call create_L0_tables before using the table driven functions."
	     << endl;
	exit(1);
      }

    if (global_prec == 53) {
      if (info > 2)
	cout << "Using double" << endl;

      L0 = approximateL0_ImaginaryNumberField_table_impl<double>(*this,n);
    }
    else {
      if (info > 2)
	cout << "Using quad_float" << endl;

      L0 = approximateL0_ImaginaryNumberField_table_impl<quad_float>(*this,n);
    }

    nterms_used[L0_REF] = n;

    if (info > 2)
      cout << "Done approxL0_imag_table - L0 = " << L0 << endl;

    return L0;
  }


  template <class TYPE> RR approximateL0_ImaginaryNumberField_table_impl(L_function<long>& lfunc, long n)
  {
    long int i,Di;	      // loop var
    TYPE T1, T2, temp1, temp2;  // used in the calulation of the two multipliers
    TYPE f, h, en;	      // used to determine en inductively
    TYPE S1, S2, S3;	      // the three sums in the calculation
    long dk,D2,D8;
    TYPE rD,CHI,CHI_n;

    /* Calculate the multipliers for the two sums:
     * T1 = sqrt(Delta) / PI
     * T2 = sqrt(PI/Delta)
     */
    dk = std::abs (lfunc.Delta);
    rD = sqrt(to<TYPE>(dk));

    T1 = rD / to<TYPE>(lfunc.PI);
    T2 = to<TYPE>(1) / rD;

    /* setup the variables for the calculation of en */
    f = exp (to<TYPE>(-lfunc.PI) / to<TYPE> (dk));
    h = f*f;
    en = f;

    clear (S1);
    clear (S2);
    clear (S3);
    clear (temp1);
    clear (temp2);

    D2 = (lfunc.Delta & 1);
    D8 = (lfunc.Delta & 7);

    for (i = 1; i < n; i++)
      {
	if (!(i & 1) && !D2) // D and p are both even
	  CHI = CHI_n = to<TYPE>(0);
	else
	  {
	    bool neg = false;
	    long t = i;
	    //          long temp;

	    while (!(t & 1))
	      {
		t >>= 1;
		neg = !neg;
	      }

	    if (t > 1 && t < lfunc.table_size) {
	      Di = lfunc.Delta % t;
	      if (Di < 0)  Di += t;
	      CHI = to<TYPE>(lfunc.chi_table[(t-1) >> 1][Di]);
	    }
	    else
	      CHI = lfunc.Chi.quadratic_long(t);

	    if (neg && (D8 == 5))
	      CHI = -CHI;

	    CHI_n = CHI/to<TYPE>(i);
	  }

	// use en in the first sum,
	temp1 = CHI_n;
	temp1 *= en;

	S1 += temp1;
	S3 += CHI;

	temp1 = en;		// S2
	// calculate en+1 for the S2 sum
	f *= h;
	en *= f;

	temp1 += en;
	temp2 = S3 * temp1;
	S2 += temp2;
      }

    lfunc.L0_result.assign (to<RR>((S1 * T1) + (S2 * T2)), to<RR>(0));

    return lfunc.L0_result.real ();
  }




  /****************************** L1 FUNCTIONS ***********************************************/

  /*
   * Function: L_function<long>::approximateL1_table
   * Description: This function calculates the value of L1 using eric bach's method
   *              of truncated euler terms. It used doubles for internal calculations
   * Inputs: long terms - The number of terms in the summation to calculate
   * Outputs: RR - The final value of the calculation.
   */
  template <> RR L_function < long >::approximateL1_table ()
  {
    double C, E, wt;
    long i,Q,Q2,P,Dp;
    long *pl;

    if (info > 2)
      cout << "\nStarting approxL1_table:  Q2 = " << (global_num_terms << 1) << ", table_size = " << table_size << endl;

    /* make sure the tables are created and ready to go */
    if (!L1table)
      {
	cerr <<
	  "Please call create_L1_tables before using the table driven functions."
	     << endl;
	exit(1);
      }

    Q = global_num_terms;
    Q2 = Q << 1;
    if (Q2 > table_size) {
      cerr << "ERROR (approximateL1_table):  global_num_terms set incorrectly" << endl;
      exit(1);
    }

    // compute weight
    clear (C);
    for (i = Q; i <= Q2 - 1; ++i)
      C += log_table[i];             // log_table[i] = i log(i)

    // compute p=2 term
    long m8 = (Delta & 7);
    if (m8 == 1)
      E = log( (double) 2.0);
    else if (m8 == 5)
      E = log( (double) 2.0 / (double) 3.0);
    else
      clear(E);

    // compute partial product  3 <= p < Q
    pl = prime_list;              // *pl = 3
    P = *pl;
    while (P < Q)
      {
	Dp = Delta % P;
	if (Dp < 0)  Dp += P;
	E += L1table[P][Dp];
	++pl;
	P = *pl;
      }

    // computed weighted partial products for Q < p < 2Q
    set (wt);
    for (i = Q; i <= P; ++i)
      wt -= log_table[i] / C;

    while (P < Q2)
      {
	Dp = Delta % P;
	if (Dp < 0)  Dp += P;
	E += wt * L1table[P][Dp];
	++pl;
	P = *pl;
	wt -= (log_table[P - 1] + log_table[P]) / C;
      }

    nterms_used[L1_REF] = global_num_terms;
    L1_result = exp(to_RR(E));

    if (info > 2)
      cout << "Done approxL1_table - L1 = " << L1_result << endl;

    return L1_result;
  }


} // ANTL
