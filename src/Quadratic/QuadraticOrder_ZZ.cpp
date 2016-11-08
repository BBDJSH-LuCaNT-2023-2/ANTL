/**
 * @file QuadraticOrder_ZZ.cpp
 * @author Michael Jacobson
 */

#include <ANTL/Quadratic/QuadraticOrder.hpp>

namespace ANTL
{
  template <>
  QuadraticOrder<ZZ>::QuadraticOrder(const ZZ & D)
  {
    // test whether D is a valid discriminant
    long m4 = rem(D,4);
    if (m4 < 0)
      m4 += 4;

    if (m4 == 0 || m4 == 1) {
      ZZ rD = SqrRoot (abs (D));
        if (rD * rD != abs (D)) {
          // assign new values
          Delta = D;
          g = 0;
          hx = 0;

	  /*
          // initialize invariant values
          clear (R);
          clear (h);
          rank = 0;
          clear (L);

          Lfunc.init(Delta,QUADRATIC_MODE);
	  */
	}
    }
  }



  //
  // QuadraticOrder<ZZ>::IsImaginary()
  //
  // Task:
  //      returns true if the quadratic order is imaginary
  //

  template <> bool QuadraticOrder<ZZ>::IsImaginary () const
  {
    return (Delta < 0);
  }



  //
  // QuadraticOrder<ZZ>::IsUnusual()
  //
  // Task:
  //      returns false, as quadratic orders are never unusual
  //

  template <> bool QuadraticOrder < ZZ >::IsUnusual () const
  {
    return false;
  }



  //
  // QuadraticOrder<ZZ>::IsReal()
  //
  // Task:
  //      returns true if the quadratic order is real
  //

  template <> bool QuadraticOrder < ZZ >::IsReal () const
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
  QuadraticOrder<ZZ> & randomImaginaryOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    long m4 = RandomWord () % 2;
    ZZ disc, rD, p;

    do
      {
	if (prime)
	  {
	    if (m4) {
	      do {
		p = RandomPrime_ZZ (size);
	      } while (((p % 4) != 1) && (NumBits(p) != size));
	      disc = p;
	    }
	    else {
	      do {
		p = RandomPrime_ZZ (size - 2);
	      } while (((p % 4) != 3) && (NumBits (4*p) != size));
	      disc = 4*p;
	    }
	  }
	else
	  disc = RandomLen_ZZ (size);
	disc = -abs (disc);
  
        rD = SqrRoot (abs (disc));
      }
    while (rD * rD == abs (disc));

    return new QuadraticOrder<ZZ>(disc);
  }



  //
  // randomUnusualOrder
  //
  // Task:
  //      generates a random real order
  //

  template <>
  QuadraticOrder<ZZ> & randomUnusualOrder (long size, bool prime)
  {
  }



  //
  // randomRealOrder
  //
  // Task:
  //      generates a random real order
  //

  template <>
  QuadraticOrder<ZZ> & randomRealOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    long m4 = RandomWord () % 2;
    ZZ disc,rd,p;

    do
      {
	if (prime)
	  {
	    if (m4) {
	      do {
		p = RandomPrime_ZZ (size);
	      } while (((p % 4) != 1) && (NumBits(p) != size));
	      disc = p;
	    }
	    else {
	      do {
		p = RandomPrime_ZZ (size - 2);
	      } while (((p % 4) != 3) && (NumBits (4*p) != size));
	      disc = 4*p;
	    }
	  }
	else
	  disc = RandomLen_ZZ (size);
	disc = abs (disc);

        rD = SqrRoot (disc);
      } 
    while (rD * rD == disc);

    return new QuadraticOrder<ZZ>(disc);
  }
*/



  /*
  //
  // QuadraticOrder<ZZ>::use_Lfunction_tables(H)
  //
  // Task:
  //      initializes table-driven L-function approximations, assuming that
  //      discriminants will be less than H
  //

  template <> void QuadraticOrder < ZZ >::use_Lfunction_tables (const ZZ & H)
  {
    use_tables = true;
    if (unconditional)
      Lfunc.create_L0_tables(H, log (sqrt (double (2))));
    else
      Lfunc.create_L1_tables(H, log (sqrt (double (2))));
  }


  //
  // QuadraticOrder<ZZ>::set_Lfunction_global(H)
  //
  // Task:
  //      sets the global terms and precision constants for L_function
  //      computations
  //
                                                                              
  template <> 
  void QuadraticOrder < ZZ >::set_Lfunction_global (const ZZ & H)
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
  // QuadraticOrder<ZZ>::Lfunction
  //
  // Task:
  //      returns an approximation of L_K(1/q) computed via the analytic class
  //      number formula.
  //

  template <> RR QuadraticOrder < ZZ >::Lfunction ()
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
  // QuadraticOrder<ZZ>::LDfunction
  //
  // Task:
  //      returns an approximation of the LD(1) function defined by Shanks
  //      (L(1,X) without the 2-factor) computed via the analytic class
  //      number formula.
  //

  template <> RR QuadraticOrder < ZZ >::LDfunction ()
  {
    RR LD;
    long m8;

    m8 = Delta % 8;
    if (m8 < 0)
      m8 += 8;
    LD = Lfunction ();
    if (m8 == 1)
      LD /= 2;
    else if (m8 == 5)
      LD *= to_RR (1.5);

    return LD;
  }



  //
  // QuadraticOrder<ZZ>::LLI
  //
  // Task:
  //      returns an approximation of the Lower Littlewood Index (defined by
  //      Shanks).
  //

  template <> RR QuadraticOrder < ZZ >::LLI ()
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

	LLI = Lfunction () * c * log (log (to_RR (abs (Delta))));
      }

    return LLI;
  }



  //
  // QuadraticOrder<ZZ>::LLI_D
  //
  // Task:
  //      returns an approximation of the LLI_D value (defined by Shanks).
  //

  template <> RR QuadraticOrder < ZZ >::LLI_D ()
  {
    RR LLI, c;

    if (::IsZero (Delta))
      clear (LLI);
    else
      {
	//    c = to_RR(8) * exp(Euler()) / (ComputePi_RR() * ComputePi_RR());
	LLI = LDfunction () * c * log (log (to_RR (abs (Delta << 2))));
      }

    return LLI;
  }



  //
  // QuadraticOrder<ZZ>::ULI
  //
  // Task:
  //      returns an approximation of the Upper Littlewood Index (defined by
  //      Shanks).
  //

  template <> RR QuadraticOrder < ZZ >::ULI ()
  {
    RR ULI, c;

    if (::IsZero (Delta))
      clear (ULI);
    else
      {
	//    c = exp(Euler());

	if (IsOdd (Delta))
	  c *= 2;

	ULI = Lfunction () / (c * log (log (to_RR (abs (Delta)))));
      }

    return ULI;
  }



  //
  // QuadraticOrder<ZZ>::ULI_D
  //
  // Task:
  //      returns an approximation of the ULI_D value (defined by Shanks).
  //

  template <> RR QuadraticOrder < ZZ >::ULI_D ()
  {
    RR ULI, c;

    if (::IsZero (Delta))
      clear (ULI);
    else
      {
	//    c = exp(Euler());
	ULI = LDfunction () / (c * log (log (to_RR (abs (Delta << 2)))));
      }

    return ULI;
  }



  //
  // QuadraticOrder<ZZ>::set_regulator(newR)
  //
  // Sets the regulator of the order - no error checking!
  //

  template <> void QuadraticOrder < ZZ >::set_regulator (const ZZ & newR)
  {
    qi_pair < ZZ > A, U;
    U.assign_one ();
    nuclose (A, newR);
    A.adjust (U);
    A.get_distance (R);
    Rbsgs = false;
    Rconditional = true;
  }


  //
  // QuadraticOrder<ZZ>::set_regulator(newR)
  //
  // Sets the regulator of the order - no error checking!
  //

  template <>
  void QuadraticOrder < ZZ >::set_regulator (const qo_distance < ZZ > &newR)
  {
    R = newR;
    Rbsgs = false;
    Rconditional = true;
  }



  // distance computation in BS-GS methods

  template <>
  void
  QuadraticOrder < ZZ >::
  combine_BSGS (qo_distance < ZZ > &dist, const qi_pair < ZZ > &DD,
		const qo_hash_entry_real < ZZ > *F)
  {
    qi_pair < ZZ > C (*F);
    C.adjust (DD);
    divide (dist, DD.get_distance (), C.get_distance ());
  }


  template <>
  void
  QuadraticOrder < ZZ >::
  combine_conj_BSGS (qo_distance < ZZ > &dist, const qi_pair < ZZ > &DD,
		     const qo_hash_entry_real < ZZ > *F)
  {
    qi_pair < ZZ > C (*F);
    C.adjust (-DD);
    multiply (dist, DD.get_distance (), C.get_distance ());
    dist.divide (DD.get_a ());
  }

  template <>
  void
  QuadraticOrder < ZZ >::
  combine_BSGS_easy (qo_distance < ZZ > &dist, const qi_pair < ZZ > &DD,
		     const qo_hash_entry_real < ZZ > *F)
  {
    qi_pair < ZZ > C (*F);
    C.adjust (DD);
    divide(dist,C.get_distance(),DD.get_distance());
  }


  template <>
  void
  QuadraticOrder < ZZ >::
  combine_conj_BSGS_easy (qo_distance < ZZ > &dist, const qi_pair < ZZ > &DD,
			  const qo_hash_entry_real < ZZ > *F)
  {
    qi_pair < ZZ > C (*F);
    C.adjust (-DD);
    divide(dist,R,C.get_distance());
    divide(dist,dist,DD.get_distance());
    dist.multiply(DD.get_a());
  }
  */


  //
  // operator <<
  //
  // Task:
  //      outputs a QuadraticOrder to the ostream out.
  //

  template <>
  void
  QuadraticOrder<ZZ>::write_to_file(std::ostream & out) const
  {
    out << "Quadratic order:" << endl;
    out << "   Delta = " << Delta << " (" << NumBits (Delta) << " bits)" << endl;

    /*
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
    */
  }


  /*
  template <>
  void QuadraticOrder<ZZ>::bsgs_getentrysize (ZZ & entry_size, bool nodist)
  {
    ZZ temp;

    entry_size = 3 * NTL_BITS_PER_LONG;	// space for pointers
    if (nodist)
      {
	entry_size += SqrRoot (Delta).size () * NTL_ZZ_NBITS;	// a coeff
	entry_size += 8;		// b coeff
      }
    else
      {
	temp = to_ZZ (1) << qo_distance < ZZ >::get_p ();
	entry_size +=
	  temp.size () * NTL_ZZ_NBITS +
	  3 * SqrRoot (Delta).size () * NTL_ZZ_NBITS;
      }

    // convert bits to bytes
    entry_size >>= 2;
  }


  template <> void QuadraticOrder < ZZ >::get_next_prime (qi_class < ZZ > &G)
  {
    static PrimeSeq PS;

    if (G.is_one ())
      PS.reset (0);

    long p = PS.next ();

    while (!G.assign_prime (to_ZZ (p)))
      p = PS.next ();
  }
  */

} // ANTL

