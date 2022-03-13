/**
 * @file QuadraticOrder_ZZ.cpp
 * @author Michael Jacobson
 */

#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <ANTL/L_function/L_function.hpp>

using namespace ANTL;

template <>
QuadraticOrder<ZZ>::QuadraticOrder(const ZZ & D)
{
  // test whether D is a valid discriminant
  long m4 = rem(D,4);
  if (m4 < 0) {
    m4 += 4;
  }

  if (m4 == 0 || m4 == 1) {
    ZZ rD = SqrRoot (abs (D));
    if (rD * rD != abs (D)) {
      // assign new values
      Delta = D;
      g = 0;
      hx = 0;

      //Setting preferred strategies. More nuanced logic should eventually go here.
      preferred_red = 0;
      preferred_mul = 0;
      preferred_sqr = 0;
      preferred_cube = 0;

      //Explicitly setting arithmetic object pointer to nullptr
      red_best = nullptr;
      red_plain_imag = nullptr;
      red_plain_real = nullptr;
      red_fast  = nullptr;

      mul_best = nullptr;
      mul_plain = nullptr;
      mul_nucomp = nullptr;

      sqr_best = nullptr;
      sqr_plain = nullptr;
      sqr_nudupl = nullptr;

      cube_best = nullptr;
      cube_plain = nullptr;
      cube_nucube = nullptr;
      cube_mulsqr = nullptr;
      /*
        // initialize invariant values
        clear (R);
        clear (h);
        rank = 0;
        clear (L);
      */
        Lfunc.init(Delta, QUADRATIC_MODE);

    }
  }
}

  //=============================
  // Beginning of qo_lfunction.cpp port
  //=============================

  //
  // QuadraticOrder<ZZ>::generate_optimal_Q()
  //
  // Task:
  //      computes the optimal value of Q for computing h* such that
  //      h* < h < 2h*
  //

  template <> long QuadraticOrder <ZZ>::generate_optimal_Q() {
    long OQ;
    RR A, l2;
    PrimeSeq primes;

    l2 = log (sqrt (to_RR (2)));

    // set OQ = 3
    OQ = primes.next ();

    do {
      OQ = primes.next ();
      A = Lfunc.calculate_L1_error (Delta, OQ);
    }
    while (A >= l2);

    return OQ;
  }

  //
  // QuadraticOrder<ZZ>::get_optimal_Q()
  //
  // Task:
  //      returns a value of Q from a pre-computed table which will compute h*
  //      such that h* < h < 2h*.
  //

  template <> long QuadraticOrder < ZZ >::get_optimal_Q() {
    long Dlog;
    ZZ temp;

    temp = FloorToZZ (log10 (to_RR (abs (Delta))));
    conv (Dlog, temp);

    if ((Dlog / 5) > 19) {
      return generate_optimal_Q ();
    }

    return OQvals[Dlog / 5];
  }

  //
  // QuadraticOrder<ZZ>::get_optimal_Q_cnum()
  //
  // Task:
  //      returns a value of Q from a pre-computed table which is sufficient
  //      to determine whether the approximation of h gleaned from the analytic
  //      class number formula is actually the class number, provided that h <=3.
  //

  template <> long QuadraticOrder < ZZ >::get_optimal_Q_cnum() {
    long Dlog;
    ZZ temp;

    temp = FloorToZZ (log10 (to_RR (abs (Delta))));
    conv (Dlog, temp);

    return OQvals_cnum[Dlog / 5];
  }

  //
  // QuadraticOrder<ZZ>::generate_optimal_Q_cfunc()
  //
  // Task:
  //      computes the optimal value of Q required to compute C(Delta) to 8
  //      significant digits.
  //

  template <> long QuadraticOrder < ZZ >::generate_optimal_Q_cfunc() {
    long QQ;
    RR OQ, B, sb, val, lDelta, lQ, temp1, temp2;

    val = log (1 + (SqrRoot (to_RR (1.0000002)) / 2));

    lDelta = log (to_RR (abs (Delta)));

    OQ = 100;

    lQ = log (OQ);
    B =
      lDelta * (inv (ComputePi_RR () * lQ) + to_RR (5.3) / (lQ * lQ)) +
      to_RR (4) / lQ;
    B += to_RR (1) / ComputePi_RR ();

    sb = 2 / (OQ * OQ) + B * (13 * lQ + 8) / (9 * OQ * SqrRoot (OQ));

    while (sb >= val) {
      if (OQ == 100) {
        OQ = 1000;
      }
      else if (OQ < 5000) {
        OQ += 1000;
      }
      else {
        OQ += 5000;
      }

      lQ = log (OQ);
      B = lDelta * (inv (ComputePi_RR () * lQ) + to_RR (5.3) / (lQ * lQ)) + to_RR (4) / lQ;
      B += to_RR (1) / ComputePi_RR ();

      sb = 2 / (OQ * OQ) + B * (13 * lQ + 8) / (9 * OQ * SqrRoot (OQ));
    }

    conv (QQ, OQ);
    return QQ;
  }

  //
  // QuadraticOrder<ZZ>::get_optimal_Q_cfunc()
  //
  // Task:
  //      returns a value of Q from a pre-computed list which is sufficient to
  //      compute an approximation of C(Delta) accurate to 8 significant digits.
  //

  template <> long QuadraticOrder < ZZ >::get_optimal_Q_cfunc() {
    long Dlog;
    ZZ temp;

    temp = FloorToZZ (log10 (to_RR (abs (Delta))));
    conv (Dlog, temp);

    if ((Dlog / 5) > 19) {
      return get_optimal_Q_cfunc ();
    }

    return OQvals_cfunc[Dlog / 5];
  }

  //
  // QuadraticOrder<ZZ>::estimate_C
  //
  // Task:
  //      computes a truncated product approximation of C(D) using primes
  //      less than Q.
  //

  template <> RR QuadraticOrder < ZZ >::estimate_C (long nQ) {
    PrimeSeq primes;
    long jac, P;
    double E;

    E = (double) 1.0;

    P = primes.next ();

    // compute partial product for LQ <= p < QQ
    while (P < nQ) {
      P = primes.next ();
      jac = Jacobi (Delta, P);
      E *= (double) 1.0 - ((double) jac / (P - 1.0));
    }

    return to_RR (E);
  }

  //
  // QuadraticOrder<ZZ>::approximate_hR()
  //
  // Task:
  //      returns an approximation of hR
  //

  template <> ZZ QuadraticOrder < ZZ >::approximate_hR() {
    ZZ hR;
    RR temp,FI;

    std::cout <<"DEBUG: approximate_hR() hello1!" << std::endl;

    long n = get_optimal_Q_cnum ();

    FI = Lfunc.approximateL1 (n);

    if (is_imaginary()) {
      // h = sqrt(Delta) L / Pi
      temp = FI * SqrRoot (to_RR (-Delta)) / ComputePi_RR ();
      if (Delta == -4) {
        temp *= 2;
      }
      if (Delta == -3) {
        temp *= 3;
      }
    }
    else {
      // hR = sqrt(Delta) L / 2
      temp = FI * SqrRoot (to_RR (Delta)) / 2;
    }

    hR = CeilToZZ (temp);
    return hR;
  }

  //
  // QuadraticOrder<ZZ>::lower_bound_hR()
  //
  // Task:
  //      returns a lower bound of hR such that L < hR < 2L
  //

  template <> ZZ QuadraticOrder < ZZ >::lower_bound_hR (){
    ZZ hR;
    RR temp,FI;

    std::cout <<"DEBUG: lower_bound_hR hello1!" << std::endl;

    if (unconditional) {
      // compute hR apporximation

      if (is_imaginary()) {
        if (use_tables) {
          FI = Lfunc.approximateL0_ImaginaryNumberField_table(log (sqrt (double (2))));
        }
        else {
          FI = Lfunc.approximateL0_ImaginaryNumberField(log (sqrt (double (2))));
        }

        temp = FI;
        if (Delta == -4)
          temp *= 2;
        if (Delta == -3)
          temp *= 3;
      }
      else {
        if (use_tables) {
          FI = Lfunc.approximateL0_RealNumberField_table(log (sqrt (double (2))));
        }
        else {
          FI = Lfunc.approximateL0_RealNumberField(log (sqrt (double (2))));
        }
        temp = FI;
      }
    }
    else {
      if (use_tables) {
        FI = Lfunc.approximateL1_table();
      }
      else {
        long n = get_optimal_Q();
        FI = Lfunc.approximateL1 (n);
      }



      if (is_imaginary())  {
        // h = sqrt(Delta) L / Pi
        temp = FI * SqrRoot (to_RR (-Delta)) / ComputePi_RR ();
        if (Delta == -4) {
          temp *= 2;
        }
        if (Delta == -3) {
          temp *= 3;
        }
      }
      else {
        // hR = sqrt(Delta) L / 2
        temp = FI * SqrRoot (to_RR (Delta)) / 2;
      }
    }

    hR = CeilToZZ (temp / SqrRoot (to_RR (2)));
    return hR;
  }



  //
  // QuadraticOrder<ZZ>::estimate_hR_error(long n)
  //
  // Task:
  //      returns L such that |hR - hR'| < exp(L)^2
  //

  template <> ZZ QuadraticOrder < ZZ >::estimate_hR_error() {
    ZZ err;
    RR Aval, Fval, temp;

    std::cout <<"DEBUG: estimate_hR_error hello1!" << std::endl;

    if (Lfunc.terms_used(1) == 0)  return ZZ::zero();

    long n = get_optimal_Q_cnum ();
    RR FI = Lfunc.approximateL1(n);

    Aval = Lfunc.calculate_L1_error (Delta, Lfunc.terms_used(1));
    Fval = exp (Aval) - 1;
    temp = 1 - exp (-Aval);
    if (temp > Fval)
      Fval = temp;

    if (is_imaginary()) {
      // h = sqrt(Delta) L / Pi
      Fval *= FI * SqrRoot (to_RR (-Delta)) / ComputePi_RR ();
      if (Delta == -4) {
        temp *= 2;
      }
      if (Delta == -3) {
        temp *= 3;
      }
    }
    else {
      // hR = sqrt(Delta) L / 2
      Fval *= FI * SqrRoot (to_RR (Delta)) / 2;
    }

    err = CeilToZZ (Fval * log (to_RR (2))) >> 1;
    return err;
  }

  //=============================
  // End of qo_lfunction.cpp port
  //=============================

  //
  // QuadraticOrder<ZZ>::is_imaginary()
  //
  // Task:
  //      returns true if the quadratic order is imaginary
  //

  template <> bool QuadraticOrder<ZZ>::is_imaginary () const
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
  // QuadraticOrder<ZZ>::is_real()
  //
  // Task:
  //      returns true if the quadratic order is real
  //

  template <> bool QuadraticOrder < ZZ >::is_real () const
  {
    return (Delta > 0);
  }

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

