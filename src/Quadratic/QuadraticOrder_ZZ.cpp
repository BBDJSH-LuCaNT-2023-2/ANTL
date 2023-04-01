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
    floor_root_delta = rD;
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
      //Lfunc.init(Delta, QUADRATIC_MODE);

    }
  }
}

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

