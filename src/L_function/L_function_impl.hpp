/**
 * @file L_function_impl.hpp
 * @author Leonard Nooy
 * @remarks This file is to be included from L_function.h only.
 * @version $Header$
 *
 * This file contains the code that the compiler will instantiate when the
 * user attempts to create a class that is not one of the types: long, ZZ,
 * GF2EX, ZZ_px, zz_px, ZZ_pEX, zz_PEX
 */

/*
 * Revision History:
 *    May 21, 2003 - File created
 *    May 28, 2003 - Added function to calculate optimal number of
 *                   terms to calculate Euler product to, based on
 *                   an error bound.
 *    June 4, 2003 - Modified this file so that it would instantiate
 *                   the code for L(1) over ZZ_pX, ZZ_pEX, zz_pX, zz_pEX.
 *  August 14 2003 - Moved the functions for ZZ_px, ZZ_pEX, etc to
 *                   L_function_FF.c file.
 *    Jan 13, 2003 - moved approximateL0_Quartic(double) here, as it is the same
 *                   in all the L_function specializations (Rennie deGraaf)
 */

#include <string.h>

namespace ANTL
{

  template < class T > long L_function < T >::info = 0;

  /*
   * Function: L_function<T>::L_function
   * Description: constructor
   * Inputs: none.
   * Outputs: none.
   */
  template < class T > L_function < T >::L_function ()
  {
    memset (nterms_used, 0, sizeof (long) * 3);
    mode = QUADRATIC_MODE;
    global_flag = false;
    global_num_terms = 0;
    global_prec = 0;

    // get the value of PI to 100 bits of accuracy
    RR::SetPrecision (100);
    ComputePi (PI);

    table_size = 0;
    log_table = 0;
    L1table = 0;
    chi_table = 0;
    prime_list = 0;
    numP = 0;

    epsilon = 0;
  }

  /*
   * Function: L_function<T>::~L_function
   * Description: Destructor
   * Inputs: none.
   * Outputs: none.
   */
  template < class T > L_function < T >::~L_function ()
  {
    delete_tables();
  }


  /*
   * Function: L_function<T>::delete_tables
   * Description: Destructor
   * Inputs: none.
   * Outputs: none.
   */
  template < class T > 
  void L_function < T >::delete_tables ()
  {
    long i;

    if (log_table)
      {
	delete [] log_table;
	log_table = 0;
      }

    if (L1table)
      {
	for (i=0; i<table_size; ++i)
	  if (L1table[i])
	    delete [] L1table[i];     
	delete [] L1table;
	L1table = 0;
      }

    if (chi_table)
      {
	for (i=0; i<table_size_real; ++i)
	  if (chi_table[i])
	    delete [] chi_table[i];     
	delete [] chi_table;
	chi_table = 0;
      }

    table_size = 0;

    if (prime_list) {
      delete [] prime_list;
      prime_list = 0;
    }

    numP = 0;
  }


  //
  // L_function::verbose()
  //
  // Task:
  //      sets the verbosity of commands.  Currently, the following levels are
  //      supported:
  //         0 - nothing
  //         1 - run times and some run-time data (data only for subexp algs.)
  //         >1 - state information for subexponential algorithms
  //

  template < class T > void L_function < T >::verbose (long state)
  {
    if (state <= 0)
      {
	info = 0;
      }
    else
      {
	info = state;
      }
  }


  /*
   * Function: L_function<T>::init
   * Description: This function must be called before using the rest of the L_function class
   * Inputs: const T & x - The conductor of the calculation to preform
   *         int m - The mode in which the class should run (This only effects the L0 calculations )
   * Outputs: none.
   */
  template < class T > void L_function < T >::init (const T & x, int m)
  {
    if (m == QUADRATIC_MODE || m == QUARTIC_MODE)
      mode = m;
    else
      {
	cout << "Init: Unrecognized Mode, Defaulting to Quadratic Mode." <<
	  endl;
	mode = QUADRATIC_MODE;
      }

    Delta.kill();
    hx.kill();
    a.kill();
    b.kill();

    Delta = x;
    memset (nterms_used, 0, sizeof (long) * 3);
    T temp = Delta;
    MakeMonic(temp);
    r2 = compute_r2(temp);

    RR::SetPrecision(53);

    if (mode == QUADRATIC_MODE)
      Chi.init(Delta,QUADRATIC_MODE);
    else
      Chi.init(Delta,QUARTIC_MODE);
  }

  /*
   * Function: L_function<T>::init
   * Description: This function fills the same purpose as the one above, but it takes two initalization
   *              varaiable instead of one because two are needed by GF2EX
   * Inputs: const T & x - The discriminant
   *         const T & y - The conductor
   *         int m - The mode to run in. (This has no meaning since calling this init function usually means
   *                 you are using the class to do calculation over Gf2EX)
   * Outputs: none.
   */
  template < class T >
  void L_function < T >::init (const T & x, const T & y, int m)
  {
    if (m == QUADRATIC_MODE || m == QUARTIC_MODE)
      mode = m;
    else
      {
	cout << "Init: Unrecognized Mode, Defaulting to Quadratic Mode." <<
	  endl;
	mode = QUADRATIC_MODE;
      }

    Delta.kill();
    hx.kill();
    a.kill();
    b.kill();

    Delta = x;
    hx = y;
    memset (nterms_used, 0, sizeof (long) * 3);

    r2 = compute_r2(Delta);

    RR::SetPrecision (53);

    if (mode == QUADRATIC_MODE)
      Chi.init(Delta,hx,QUADRATIC_MODE);
    else
      Chi.init(Delta,hx,QUARTIC_MODE);
  }

  /*
   * Function::L_function<T>::set_global_prec
   * Description: This function takes an upper bound on a discriminant range and the runtime mode
   *              and calculates an appropriate precision to ensure the round-off error will not
   *              become a problem.
   * Inputs: const T & D - The discriminant upper bound
   *         int m - The mode (QUARTIC_MODE or QUADRATIC_MODE)
   * Outputs: none.
   */
  template < class T >
  void L_function < T >::set_global_prec (const T & D, int m)
  {
    if (m == QUADRATIC_L0_MODE)
      {
	global_flag = true;
	if (D > 0)
	  {			// real quadratic number field
	    global_prec = calculate_precision_RealNumberField (D);
	    RR::SetPrecision (global_prec);
	  }
	else
	  {			// imaginary quadratic number field
	    global_prec = calculate_precision_ImaginaryNumberField (D);
	    RR::SetPrecision (global_prec);
	  }
      }
    else if (m == QUADRATIC_L1_MODE)
      {
	global_flag = true;
	global_prec = calculate_precision (D);
	RR::SetPrecision (global_prec);
      }
    else if (m == QUARTIC_MODE)
      {
	global_flag = true;
	global_prec = calculate_precision_Quartic (D);
	RR::SetPrecision (global_prec);

	/* This function is usually called when calculating large tables of h- numbers
	   so pass the upper bound on the discriminant to the quartic character stuff
	   that way it can allocate one big array that doesn't get deallocated
	*/
	Chi.allocate_quartic_table();
      }
    else
      cout <<
	"QUARTIC_MODE, QUADRATIC_L1_MODE and QUADRATIC_L0_MODE are the only modes supported at this time\n";
  }

  /*
   * Function:: L_function<T>::set_global_terms
   * Description: This function takes an upper bound on a discriminant range and the runtime mode
   *              and calculates an appropriate number of terms to ensure the round-off errror will
   *              not become a problem.
   * Inputs: const T & D - The discriminant upper bound
   *         double acc - An upper bound on the absolute error allowable by the calculation
   *         int m - The mode (QUARTIC_MODE or QUADRATIC_L0_MODE or QUADRATIC_L1_MODE)
   * Outputs: none.
   */
  template < class T >
  void L_function < T >::set_global_terms (const T & D, double acc, int m)
  {
    if (m == QUADRATIC_L0_MODE)
      {
	global_flag = true;
	if (D > 0)
	  {			// real number field
	    global_num_terms = calculate_optimal_n_RealNumberField (D, acc);
	  }
	else
	  {			// imaginary number field
	    global_num_terms =
	      calculate_optimal_n_ImaginaryNumberField (D, acc);
	  }
      }
    else if (m == QUADRATIC_L1_MODE)
      {
	global_flag = true;
	global_num_terms = calculate_optimal_terms_new (D, acc);
      }
    else if (m == QUARTIC_MODE)
      {				// imaginary cyclic quartic number field
	global_num_terms = calculate_optimal_n_Quartic (D, acc);
      }
    else
      cout <<
	"QUARTIC_MODE,QUADRATIC_L1_MODE and QUADRATIC_L0_MODE are the only modes supported at this time\n";
  }

  /*****************************************************************************************

MISC FUNCTIONS

  *****************************************************************************************/

  /*
   * Function: L_function<long>::approximateL0_Quartic
   * Description: This function is simply a wrapper for the calculate_optimal_n and approximateL0(long n)
   * Inputs: double acc - The required accuracy
   * Outputs: CC - The result of the calculation
   */
  template<class T> CC<RR> L_function<T>::approximateL0_Quartic_default(double acc)
  {
    long Q;
    long prec;

#ifndef USE_CHI2
    Chi.init_quartic_table ();
#endif

    if (global_flag)
      {
        Q = global_num_terms;
        prec = global_prec;

        if (info > 1)
	  cout << "prec: " << prec << " Q: " << Q << endl;

#ifdef FORCE_FLOAT
        approximateL0_Quartic_impl<float>(*this, Q);
#elif defined FORCE_DOUBLE
        approximateL0_Quartic_impl<double>(*this, Q);
#elif defined FORCE_XDOUBLE
        approximateL0_Quartic_impl<xdouble>(*this, Q);
#elif defined FORCE_QUAD
        approximateL0_Quartic_impl<quad_float>(*this, Q);
#elif defined FORCE_RR
        approximateL0_Quartic_impl<RR>(*this, Q);
#else
        if (prec <= 53)
	  {
            approximateL0_Quartic_impl<double>(*this, Q);
	  }
        else if (prec <= 104)
	  {
            approximateL0_Quartic_impl<quad_float>(*this, Q);
	  }
        else
	  {
            /* there is no need to set the RR precision because it was set when the global
               precision was setup */
            approximateL0_Quartic_impl<RR>(*this, Q);
	  }
#endif
      }
    else
      {
        Q = calculate_optimal_n_Quartic (Delta, acc);
        if (Q <= nterms_used[L0_REF])
	  {
            /* if the L0 function has been calculated before, return that value */
            return L0_result;
	  }

        prec = calculate_precision_Quartic (Delta);

        if (info > 1)
	  cout << "prec: " << prec << " Q: " << Q << endl;

#ifdef FORCE_FLOAT
        approximateL0_Quartic_impl<float>(*this, Q);
#elif defined FORCE_DOUBLE
        approximateL0_Quartic_impl<double>(*this, Q);
#elif defined FORCE_XDOUBLE
        approximateL0_Quartic_impl<xdouble>(*this, Q);
#elif defined FORCE_QUAD
        approximateL0_Quartic_impl<quad_float>(*this, Q);
#elif defined FORCE_RR
        RR::SetPrecision(prec);
        approximateL0_Quartic_impl<RR>(*this, Q);
#else
        if (prec <= 53)
	  {
            approximateL0_Quartic_impl<double>(*this, Q);
	  }
        else if (prec <= 104)
	  {
            approximateL0_Quartic_impl<quad_float>(*this, Q);
	  }
        else
	  {
            RR::SetPrecision (prec);
            approximateL0_Quartic_impl<RR>(*this, Q);
	  }
#endif
      }

    return L0_result;
  }

  /*
   * Function: L_function<ZZ>::approximateL0_ImaginaryNumberField
   * Description: This function is simply a wrapper for the calculate_optimal_n and approximateL0(long n)
   * Inputs: double acc - The required accuracy
   * Outputs: CC - The result of the calculation
   */
  template<class T> RR L_function<T>::approximateL0_ImaginaryNumberField_default(double acc)
  {
    long Q;
    long prec;

    if (global_flag)
      {
        Q = global_num_terms;
        prec = global_prec;

        if (prec <= 53)
	  {
            approximateL0_ImaginaryNumberField_impl<double>(*this, Q);
	  }
        else if (prec <= 104)
	  {
            approximateL0_ImaginaryNumberField_impl<quad_float>(*this, Q);
	  }
        else
	  {
            /* there is no need to set the RR precision because it was set when the global
               precision was setup */
            approximateL0_ImaginaryNumberField_impl<RR>(*this, Q);
	  }

      }
    else
      {
        Q = calculate_optimal_n_ImaginaryNumberField (Delta, acc);
        if (Q <= nterms_used[L0_REF])
	  {			/* if the L0 function has been calculated before, return that value */
            return L0_result.real ();
	  }

        prec = calculate_precision_ImaginaryNumberField (Delta);
        if (prec <= 53)
	  {
            approximateL0_ImaginaryNumberField_impl<double>(*this, Q);
	  }
        else if (prec <= 104)
	  {
            approximateL0_ImaginaryNumberField_impl<quad_float>(*this, Q);
	  }
        else
	  {
            RR::SetPrecision (prec);
            approximateL0_ImaginaryNumberField_impl<RR>(*this, Q);
	  }
      }

    return L0_result.real ();
  }

  /*
   * Function: L_function<long>::approximateL0_RealNumberField
   * Description: This function is simply a wrapper for the calculate_optimal_n and approximateL0(long n)
   * Inputs: double acc - The required accuracy
   * Outputs: CC - The result of the calculation
   */
  template<class T> RR L_function<T>::approximateL0_RealNumberField_default(double acc)
  {
    long Q;
    long prec;

    if (global_flag)
      {
        Q = global_num_terms;
        prec = global_prec;

        if (prec <= 53)
	  {
            approximateL0_RealNumberField_impl<double>(*this, Q);
	  }
        else if (prec <= 104)
	  {
            approximateL0_RealNumberField_impl<quad_float>(*this, Q);
	  }
        else
	  {
            /* there is no need to set the RR precision because it was set when the global
               precision was setup */
            approximateL0_RealNumberField_impl<RR>(*this, Q);
	  }

      }
    else
      {
        Q = calculate_optimal_n_RealNumberField (Delta, acc);
        if (Q <= nterms_used[L0_REF])
	  {			/* if the L0 function has been calculated before, return that value */
            return L0_result.real ();
	  }

        prec = calculate_precision_RealNumberField (Delta);
        if (prec <= 53)
	  {
            approximateL0_RealNumberField_impl<double>(*this, Q);
	  }
        else if (prec <= 104)
	  {
            approximateL0_RealNumberField_impl<quad_float>(*this, Q);
	  }
        else
	  {
            RR::SetPrecision (prec);
            approximateL0_RealNumberField_impl<RR>(*this, Q);
	  }
      }

    return L0_result.real ();
  }

  /* for the time being this is function only calculates using the quad_float base type */
  template<class T> CC<RR> L_function<T>::approximateL0_default(long terms)
  {
    long prec;

    if (global_prec)
      {
        if (mode == QUADRATIC_MODE)
	  {
	    if (Delta > 0)
	      {
                approximateL0_RealNumberField (0.0);
	      }
            else
	      {
                approximateL0_ImaginaryNumberField (0.0);
	      }
	  }
        else
	  {
            approximateL0_Quartic (0.0);
	  }
      }
    else
      {
        if (mode == QUADRATIC_MODE)
	  {
            if (Delta > 0)
	      {
                prec = calculate_precision_RealNumberField (Delta);
                if (prec <= 53)
		  {
                    approximateL0_RealNumberField_impl<double>(*this, terms);
		  }
                else if (prec <= 104)
		  {
                    approximateL0_RealNumberField_impl<quad_float>(*this, terms);
		  }
                else
		  {
                    RR::SetPrecision (prec);
                    approximateL0_RealNumberField_impl<RR>(*this, terms);
		  }
	      }
            else
	      {
                prec = calculate_precision_ImaginaryNumberField (Delta);
                if (prec <= 53)
		  {
                    approximateL0_ImaginaryNumberField_impl<double>(*this, terms);
		  }
                else if (prec <= 104)
		  {
                    approximateL0_ImaginaryNumberField_impl<quad_float>(*this, terms);
		  }
                else
		  {
                    RR::SetPrecision (prec);
                    approximateL0_ImaginaryNumberField_impl<RR>(*this, terms);
		  }
	      }
	  }
        else
	  {
            prec = calculate_precision_Quartic (Delta);
            if (prec <= 53)
	      {
                approximateL0_Quartic_impl<double>(*this, terms);
	      }
            else if (prec <= 104)
	      {
                approximateL0_Quartic_impl<quad_float>(*this, terms);
	      }
            else
	      {
                RR::SetPrecision (prec);
                approximateL0_Quartic_impl<RR>(*this, terms);
	      }
	  }
      }

    return L0_result;
  }

  template<class T> CC<RR> L_function<T>::approximateL0_default(double acc)
  {
    if (mode == QUADRATIC_MODE)
      {
        if (Delta > 0)
	  {
            approximateL0_RealNumberField (acc);
	  }
        else
	  {
            approximateL0_ImaginaryNumberField (acc);
	  }
      }
    else if (mode == QUARTIC_MODE)
      {
        approximateL0_Quartic (acc);
      }

    return L0_result;
  }




  /* return the first factor of the preset quartic discriminant */
  template < class T > T L_function < T >::get_a (void) const
  {
    return a;
  }

  /* return the second factor of the preset quartic discriminant */
  template < class T > T L_function < T >::get_b (void) const
  {
    return b;
  }

  /*
   * Function: L_function<T>::calculate_optimal_terms
   * Description: This function is used to figure out how many terms in
   *              the truncated Euler product will result in an error less than
   *              or equal to the given double acc.
   * Inputs: double acc - The required accuracy of the calculation
   * Outputs: long - The number of terms to use
   */
  template < class T >
  long L_function < T >::calculate_optimal_terms (const T & D, double acc)
  {
    long Q, i;
    double A, B, lq, ln, C;

    primes.reset(2);
    Q = primes.next ();
    do
      {
	i = 0;
	while (i < 16)
	  {
	    Q = primes.next ();
	    i++;
	  }
	Bach_Table (Q, A, B);

	/* calculate C */
	lq = std::log ((double) Q);
	ln = to_double (log (to<RR> (abs (D))));

	C = (A * ln + B) / (lq * std::sqrt ((double) Q));
      }
    while (acc <= C);

    return Q;
  }


  /*
   * Function: L_function<T>::calculate_L1_error
   * Description: This function is used to figure out how many terms in
   *              the truncated Euler product will result in an error less than
   *              or equal to the given double acc.
   * Inputs: double acc - The required accuracy of the calculation
   * Outputs: long - The number of terms to use
   */
  template < class T >
  RR L_function < T >::calculate_L1_error (const T & D, long Q)
  {
    double dQ = double(Q);
    double lQ,U,G,H,A,c,l;

    c = eval_c(to<RR>(abs(D)));
    l = (double(2) * (std::sqrt(double(8)) - 1)) / double(3);

    lQ = std::log(dQ);
    U = eval_U(dQ);
    G = eval_G(dQ,lQ,U,l);
    H = eval_H(dQ,lQ,U,l);

    A = c*G + H;

    return to_RR(A);
  }



  /*
   * Function: L_function<T>::calculate_optimal_terms_new
   * Description: This function is used to figure out how many terms in
   *              the truncated Euler product will result in an error less than
   *              or equal to the given double acc.
   * Inputs: double acc - The required accuracy of the calculation
   * Outputs: long - The number of terms to use
   */
  template < class T >
  long L_function < T >::calculate_optimal_terms_new (const T & D, double acc)
  {
    long Q, i;
    double A;

    primes.reset(2);
    Q = primes.next ();
    do {
      i = 0;
      while (i < 16) {
        Q = primes.next ();
        i++;
      }

    conv(A,calculate_L1_error(D,Q));
    }
    while (acc <= A);

    return Q;
  }



  template < class T > long L_function < T >::calculate_precision (const T & D)
  {
    return 53; // double precision is sufficient for almost but the largest discriminants, ie 10^1000, etc
  }

  /*
   * Function: L_function<T>::calculate_optimal_n_Quartic
   * Description: This function calculates how many terms are needed in the L0 summation
   *              to obtain an error less than or equal to the given accuracy
   * Inputs: ZZ D - The delta to calculate how many terms to use for.
   *         double acc - The required accuracy
   * Outputs: long - The number of terms to use
   */
  template < class T >
  long L_function < T >::calculate_optimal_n_Quartic (const T & D, double acc)
  {
    RR M;
    RR N;

    M = to_RR(3) / (to_RR(2)*(to_RR(acc) - to_RR(3)/(to_RR(8)*PI)));
    M = M*M / PI;

    // N > B(1/2, M, Delta) = sqrt(1/2*Delta/PI * log(M * Depta /Pi))
    N = to<RR>(D) / (to_RR(2)*PI);
    N *= log (M * to<RR>(D) / PI);
    N = ceil(SqrRoot (N));

    return to_long (N);
  }

  template < class T >
  long L_function < T >::calculate_precision_Quartic (const T & D)
  {
    /* use double precision for everything below 10^9 and quad_float for everything above */
    if (D < to_ZZ ("1000000000"))
      return 53;
    else
      return 104;
  }

  /*
   * Function: L_function<T>::calculate_optimal_n_RealNumberField
   * Description: This function calculates how many terms are needed in the L0 summation
   *              to obtain an error less than or equal to the given accuracy
   * Inputs: ZZ D - The delta to calculate how many terms to use for.
   *         double acc - The required accuracy
   * Outputs: long - The number of terms to use
   */
  template < class T >
  long
  L_function < T >::calculate_optimal_n_RealNumberField (const T & D, double acc)
  {
    RR M;
    RR N;
    RR dk;
    RR alpha,temp;


    dk = to<RR>(D);
    //  alpha = sqrt(PI/dk);
    RR rD = sqrt(dk);
    RR rPI = sqrt(PI);
    // temp = log(exp(to_RR(1))*exp(to_RR(1))/alpha);
    temp = log(exp(to_RR(1))*exp(to_RR(1))*rD/rPI);

    M = log(to_RR(2)*temp / (to_RR(4)*to_RR(acc) - (rPI*temp)/rD));

    // N > B(1/2, M, |Delta|) = sqrt((|Delta|*(log(Delta/PI) + 2M) / 2PI)
    N = dk * (log (dk / PI) + to_RR(2) * M);
    N /= (to_RR(2) * PI);
    N = ceil(SqrRoot (N));

    return to_long (N);
  }

  template < class T >
  long L_function < T >::calculate_precision_RealNumberField (const T & D)
  {
    /* use double precision for everything below 10^8 and quad_float for everything above */
    if (D < to_ZZ ("100000000"))
      return 53;
    else
      return 104;
  }

  /*
   * Function: L_function<T>::calculate_optimal_n_ImaginaryNumberField
   * Description: This function calculates how many terms are needed in the L0 summation
   *              to obtain an error less than or equal to the given accuracy
   * Inputs: ZZ D - The delta to calculate how many terms to use for.
   *         double acc - The required accuracy
   * Outputs: long - The number of terms to use
   */
  template < class T >
  long
  L_function < T >::calculate_optimal_n_ImaginaryNumberField (const T & D,
							      double acc)
  {
    RR M;
    RR N;
    RR dk;

    dk = to<RR> (abs (D));

    M = to_RR(12)*SqrRoot(dk / PI);
    M /= (to_RR(8)*to_RR(acc)*SqrRoot(dk) - to_RR(3));
    M = log(M);

    // N > B(1/2, M, |Delta|) = sqrt((|Delta|*(log(Delta/PI) + 2M) / 2PI)
    N = dk * (log (dk / PI) + to_RR(2) * M);
    N /= (to_RR(2) * PI);
    N = ceil(SqrRoot (N));

    return to_long (N);
  }

  template < class T >
  long L_function < T >::calculate_precision_ImaginaryNumberField (const T & D)
  {
    /* use double precision for everything below 10^8 and quad_float for everything above */
    if (D < to_ZZ ("100000000"))
      return 53;
    else
      return 104;
  }

  template < class T > long L_function < T >::get_epsilon (void) const
  {
    return epsilon;
  }


  /*
   * Function: L_function<ZZ>::create_L1_tables
   * Description: This function must be called before using any of the L0 table driven functions
   *              It creates and sets up the tables used by the table driven functions.
   * Inputs: ZZ D - An upper bound on the Delta to calculate
   *         double acc - The absolute error of each of the L0 functions calculations.
   * Outputs: void
   */
  template < class T >
  void L_function < T >::create_L1_tables (const T & D, double acc)
  {
    long P,Q2,a,a2;
    double val;

    if (info > 2) 
      cout << "\nStarting create_L1_tables" << endl;

    // set global precision
    set_global_prec(D,QUADRATIC_L1_MODE);

    // compute number of terms to use, assuming D is an upper bound
    set_global_terms(D,acc,QUADRATIC_L1_MODE);
    Q2 = global_num_terms << 1;

    if (info > 2)
      cout << "Q = " << global_num_terms << ", Q2 = " << Q2 << endl;

    if (Q2 > table_size) {
      // allocate tables
      delete_tables();

      table_size_real = table_size = Q2;

      L1table = new double *[Q2];
      if (!L1table) {
	cerr << "create_L1_tables - Unable to allocate memory for L1 table" << endl;
	exit(1);
      }  
      memset (L1table, 0, sizeof (double *) * Q2);

      log_table = new double[Q2];
      if (!log_table)
	{
	  cerr << "create_L1_tables - Unable to allocate memory for log tables" << endl;
	  exit(1);
	}
      memset (log_table, 0, sizeof (double) * Q2);

      prime_list = new long[Q2];
      if (!prime_list)
	{
	  cerr << "create_L1_tables - Unable to allocate memory for prime list" << endl;
	  exit(1);
	}
      memset (prime_list, 0, sizeof (long) * Q2);

      // calculate L1table
      numP = 0;
      primes.reset(3);
      P = primes.next();     // P = 3
      while (P < Q2)
	{
	  L1table[P] = new double[P];
	  if (!L1table[P])
	    {
	      cerr << "create_L1_tables - Unable to allocate memory for L1table[" << P << "]" << endl;
	      exit(1);
	    }

	  prime_list[numP] = P;
	  ++numP;

	  val = std::log((double) P / ( (double) P-1));
	  L1table[P][0] = (double) 0;
	  L1table[P][1] = val;
	  for (a=2; a<P; ++a)
	    L1table[P][a] = std::log((double) P / ((double) P + 1));

	  for (a=2; a<P; ++a) {
	    a2 = MulMod(a,a,P);
	    L1table[P][a2] = val;
	  }

	  P = primes.next ();
	}

      prime_list[numP] = P;
      ++numP;

      // calculate log table
      for (a=1; a<Q2; ++a)
	log_table[a] = std::log ((double) a) * ((double) a);
    }

    if (info > 2)
      cout << "Done create_L1_tables, table_size = " << table_size << endl;
  }



  template < class T >
  void L_function < T >::create_L0_tables(const T & D, double acc)
  {
    long i,a;

    if (info > 2) 
      cout << "\nStarting create_L0_tables" << endl;

    // set global precision
    set_global_prec(D,QUADRATIC_L0_MODE);

    // compute number of terms to use, assuming D is an upper bound
    set_global_terms(D,acc,QUADRATIC_L0_MODE);

    if (info > 2)
      cout << "n = " << global_num_terms << endl;

    if (global_num_terms > table_size && table_size < MAX_L0_TABLE) {
      // allocate tables
      delete_tables();

      table_size = global_num_terms;
      if (table_size > MAX_L0_TABLE)
	table_size = MAX_L0_TABLE;

      table_size_real = (table_size + 1) >> 1;
      chi_table = new short *[table_size_real];
      if (!chi_table) {
	cerr << "create_L0_tables - Unable to allocate memory for chi table" << endl;
	exit(1);
      }  
      memset (chi_table, 0, sizeof (short *) * (table_size_real));

      // calculate chi_table
      chi_table[1] = new short[1];
      chi_table[1][0] = (short) 1;

      for (i=3; i<table_size; i+=2) {
	chi_table[(i-1) >> 1] = new short[i];
	if (!chi_table[(i-1) >> 1])
	  {
	    cerr << "create_L0_tables - Unable to allocate memory for chi_table[" << i << "]" << endl;
	    exit(1);
	  }

	chi_table[(i-1) >> 1][0] = (short) 0;
	chi_table[(i-1) >> 1][1] = (short) 1;

	for (a=2; a<i; ++a)
	  chi_table[(i-1) >> 1][a] = (short) Kronecker(a,i);
      }
    }

    if (info > 2)
      cout << "Done create_L0_tables" << endl;
  }



  /**********************************************************************************

Functions which have no default meaning

  ***********************************************************************************/

  /*
    template < class T >
    void L_function < T >::create_L1_tables (const T & D, double acc)
    {
    cout << "create_L1_tables: Unknown Type\n";
    }

    template < class T >
    void L_function < T >::create_L0_tables (const T & D, double acc)
    {
    cout << "create_L0_tables: Unknown Type\n";
    }
  */

  template < class T > CC < RR > L_function < T >::approximateL0 (long n)
  {
    cout << "approximateL0: Unknown Type\n";
    CC < RR > a;
    a.clear ();
    return a;
  }

  template < class T >
  CC < RR > L_function < T >::approximateL0_Quartic (double acc)
  {
    cout << "approximateL0_Quartic: Unknown Type\n";
    CC < RR > a;
    a.clear ();
    return a;
  }

  template < class T >
  RR L_function < T >::approximateL0_RealNumberField (double acc)
  {
    cout << "approximateL0_RealNumberField: Unknown Type\n";
    return to_RR (0.0);
  }

  template < class T >
  RR L_function < T >::approximateL0_ImaginaryNumberField (double acc)
  {
    cout << "approximateL0_ImaginaryNumberField: Unknown Type\n";
    return to_RR (0.0);
  }

} // ANTL
