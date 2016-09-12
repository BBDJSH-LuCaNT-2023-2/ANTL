/**
 * @file QuadraticOrder_long.cpp
 * @author Michael Jacobson
 */

#include <ANTL/Quadratic/QuadraticOrder.hpp>

namespace ANTL
{
  template <>
  QuadraticOrder<long>::QuadraticOrder(const long & D)
  {
    // test whether D is a valid discriminant
    long m4 = D & 3;
    if (m4 == 0 || m4 == 1) {
      long rD = SqrRoot (labs (D));
        if (rD * rD != labs (D)) {
          // assign new values
          Delta = D;
          g = 0;
          hx = 0;

	  /*
          // initialize invariant values
          clear (R);
          clear (h);
          CL.clear ();
          gens.clear();
          clear (L);

	  Lfunc.init(Delta,QUADRATIC_MODE);
	  */
	}
    }
  }




  //
  // QuadraticOrder<long>::IsImaginary()
  //
  // Task:
  //      returns true if the function field is imaginary
  //

  template <> bool QuadraticOrder < long >::IsImaginary () const
  {
    return (Delta < 0);
  }



  //
  // QuadraticOrder<long>::IsUnusual()
  //
  // Task:
  //      returns false, as quadratic number fields are never unusual //

  template <> bool QuadraticOrder < long >::IsUnusual () const
  {
    return false;
  }




  //
  // QuadraticOrder<long>::IsReal()
  //
  // Task:
  //      returns true if the function field is real
  //

  template <> bool QuadraticOrder < long >::IsReal () const
  {
    return (Delta > 0);
  }



  /*
  //
  // randomImaginaryOrder
  //
  // Task:
  //      generates a random imaginary order
  //

  template <>
  QuadraticOrder<long> & randomImaginaryOrder (long size, bool prime)
  {
    // select a random discriminant with given genus
    long m4 = RandomWord () & 1;
    long disc,rD;
    ZZ zdisc;

    do
      {
	if (prime)
	  {
	    if (m4)
	      zdisc = RandomPrime_ZZ (size);
	    else
	      do
		zdisc = 4 * RandomPrime_ZZ (size - 2);
	      while (NumBits (zdisc) != size);
	  }
	else
	  zdisc = RandomLen_ZZ (size);
	conv (disc, -abs (zdisc));

        rD = SqrRoot (labs (disc));
      }
    while (rD * rD == labs (disc));

    return new QuadraticOrder<long>(disc);
  }


  //
  // randomUnusualOrder
  //
  // Task:
  //      generates a random real order
  //

  template <>
  QuadraticOrder<long> & randomUnusualOrder (long size, bool prime)
  {
  }



  //
  // randomRealOrder
  //
  // Task:
  //      generates a random real order
  //

  template <>
  QuadraticOrder<long> & randomRealOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    long m4 = RandomWord () & 1;
    long disc,rD;
    ZZ zdisc;

    do
      {
	if (prime)
	  {
	    if (m4)
	      zdisc = RandomPrime_ZZ (size);
	    else
	      do
		zdisc = 4 * RandomPrime_ZZ (size - 2);
	      while (NumBits (zdisc) != size);
	  }
	else
	  zdisc = RandomLen_ZZ (size);
	conv (disc, abs (zdisc));

        rD = SqrRoot (labs (disc));
      }
    while (rD * rD == labs (disc));

    return new QuadraticOrder<long>(disc);
  }
  */


  /*
  //
  // QuadraticOrder<long>::use_Lfunction_tables(H)
  //
  // Task:
  //      initializes table-driven L-function approximations, assuming that
  //      discriminants will be less than H
  //
                                                                              
  template <> void QuadraticOrder < long >::use_Lfunction_tables (const long & H)
  {
    use_tables = true;
    if (unconditional)
      Lfunc.create_L0_tables(H, log (sqrt (double (2))));
    else
      Lfunc.create_L1_tables(H, log (sqrt (double (2))));
  }



  //
  // QuadraticOrder<long>::set_Lfunction_global(H)
  //
  // Task:
  //      sets the global terms and precision constants for L_function
  //      computations
  //
                                                                              
  template <> 
  void QuadraticOrder < long >::set_Lfunction_global (const long & H)
  {
    double acc = std::log (std::sqrt (double (2)));
      
    if (unconditional) {
      // set global precision
      Lfunc.set_global_prec(H,QUADRATIC_L0_MODE);

      // compute number of terms to use, assuming D is an upper bound
      Lfunc.set_global_terms(H,acc,QUADRATIC_L0_MODE);
    }
    else {
      // set global precision
      Lfunc.set_global_prec(H,QUADRATIC_L1_MODE);

      // compute number of terms to use, assuming D is an upper bound
      Lfunc.set_global_terms(H,acc,QUADRATIC_L1_MODE);
    }
  }
                                                                              
      

  //
  // QuadraticOrder<long>::Lfunction
  //
  // Task:
  //      returns an approximation of L_K(1/q) computed via the analytic class
  //      number formula.
  //

  template <> RR QuadraticOrder < long >::Lfunction ()
  {
    if (!is_L_computed ())
      {
	// L has not been computed - compute it

	if (is_imaginary ())
	  {
	    // L = pi*h / sqrt(Delta)
	    L =
	      ComputePi_RR () * to_RR (class_number ()) /
	      SqrRoot (to_RR (-Delta));
	    if (Delta == -4)
	      L /= 2;
	    if (Delta == -3)
	      L /= 3;
	  }
	else
	  {
	    // L = 2hR / sqrt(Delta)
	    L =
	      to_RR (class_number () << 1) * regulator ().get_log () /
	      SqrRoot (to_RR (Delta));
	  }
      }

    return L;
  }



  //
  // QuadraticOrder<long>::LDfunction
  //
  // Task:
  //      returns an approximation of the LD(1) function defined by Shanks
  //      (L(1,X) without the 2-factor) computed via the analytic class
  //      number formula.
  //

  template <> RR QuadraticOrder < long >::LDfunction ()
  {
    RR LD;
    long m8;

    m8 = Delta & 7;
    LD = Lfunction ();
    if (m8 == 1)
      LD /= 2;
    else if (m8 == 5)
      LD *= to_RR (1.5);

    return LD;
  }



  //
  // QuadraticOrder<long>::LLI
  //
  // Task:
  //      returns an approximation of the Lower Littlewood Index (defined by
  //      Shanks).
  //

  template <> RR QuadraticOrder < long >::LLI ()
  {
    RR LLI, c;

    if (::IsZero (Delta))
      clear (LLI);
    else
      {
	//    c = exp(Euler()) / (ComputePi_RR() * ComputePi_RR());

	if (IsOdd (Delta))
	  c *= to_RR (12);
	else
	  c *= to_RR (8);

	LLI = Lfunction () * c * log (log (to_RR (labs (Delta))));
      }

    return LLI;
  }



  //
  // QuadraticOrder<long>::LLI_D
  //
  // Task:
  //      returns an approximation of the LLI_D value (defined by Shanks).
  //

  template <> RR QuadraticOrder < long >::LLI_D ()
  {
    RR LLI, c;

    if (::IsZero (Delta))
      clear (LLI);
    else
      {
	//    c = to_RR(8) * exp(Euler()) / (ComputePi_RR() * ComputePi_RR());
	LLI = LDfunction () * c * log (log (to_RR (labs (Delta << 2))));
      }

    return LLI;
  }



  //
  // QuadraticOrder<long>::ULI
  //
  // Task:
  //      returns an approximation of the Upper Littlewood Index (defined by
  //      Shanks).
  //

  template <> RR QuadraticOrder < long >::ULI ()
  {
    RR ULI, c;

    if (::IsZero (Delta))
      clear (ULI);
    else
      {
	//    c = exp(Euler());

	if (IsOdd (Delta))
	  c *= 2;

	ULI = Lfunction () / (c * log (log (to_RR (labs (Delta)))));
      }

    return ULI;
  }



  //
  // QuadraticOrder<long>::ULI_D
  //
  // Task:
  //      returns an approximation of the ULI_D value (defined by Shanks).
  //

  template <> RR QuadraticOrder < long >::ULI_D ()
  {
    RR ULI, c;

    if (::IsZero (Delta))
      clear (ULI);
    else
      {
	//    c = exp(Euler());
	ULI = LDfunction () / (c * log (log (to_RR (labs (Delta << 2)))));
      }

    return ULI;
  }



  //
  // QuadraticOrder<long>::set_regulator(newR)
  //
  // Sets the regulator of the order - no error checking!
  //

  template <> void QuadraticOrder < long >::set_regulator (const ZZ & newR)
  {
    qi_pair < long >A, U;

    U.assign_one ();
    nuclose (A, newR);
    A.adjust (U);
    A.get_distance (R);
    Rbsgs = false;
    Rconditional = true;
  }


  //
  // QuadraticOrder<long>::set_regulator(newR)
  //
  // Sets the regulator of the order - no error checking!
  //

  template <>
  void
  QuadraticOrder < long >::set_regulator (const qo_distance < long >&newR)
  {
    R = newR;
    Rbsgs = false;
    Rconditional = true;
  }
  */


  //
  // operator <<
  //
  // Task:
  //      outputs a QuadraticOrder to the ostream out.
  //

  /*
  template <>
  void
  QuadraticOrder<long>::write_to_file(std::ostream & out) const
  {
    out << "Quadratic order:" << endl;
    out << "   Delta = " << Delta << " (" << NumBits(Delta) << " bits)" << endl;

    if (!is_zero ())
      {
	if (is_R_computed () && is_real ())
	  out << "   R = " << R << endl;
	if (is_h_computed ())
	  out << "   h = " << h << endl;
	if (is_CL_computed ())
	  out << "   CL = " << CL << endl;
	if (rank > 0)
	  {
	    out << "   generators:" << endl;
	    for (long i = 0; i < rank; ++i)
	      {
		out << "      " << gens[i];
		if (i < rank - 1)
		  out << ",";
		out << endl;
	      }
	  }
	if (is_L_computed ())
	  out << "   L(1) = " << L << endl;
      }
  }
*/

} // ANTL

