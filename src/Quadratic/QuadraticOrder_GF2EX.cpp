/**
 * @file QuadraticOrder_GF2EX.cpp
 * @author Michael Jacobson
 */

#include <ANTL/Quadratic/QuadraticOrder.hpp>

namespace ANTL
{
  //
  // QuadraticOrder<GF2EX>::assign(const T & f, const T & h)
  //
  // Task:
  //    set to the quadratic order defined by f and h.
  //

  template <> 
  QuadraticOrder<GF2EX>::QuadraticOrder(const GF2EX & newf)
  {
  }


  template <>
  QuadraticOrder<GF2EX>::QuadraticOrder(const GF2EX & newf, const GF2EX & newh)
  {
    // degree must be >= 3
    if (deg (newf) >= 3) {
      // assign new values
      Delta = newf;
      hx = newh;

      // compute genus
      if (IsReal()){
        g = deg(hx) - 1;
      }
      else {
        if (deg(Delta) & 1){
  	  g = (deg(Delta) - 1) >> 1;
        }
        else{
	  //MDV      g = deg(hx) - 1;
	  g = (deg(Delta) >> 1) -1;
        }
      }

      /*
      // initialize invariant values
      clear (R);
      clear (h);
      gens = (qi_class < GF2EX > *)NULL;
      rank = 0;
      clear (L);

      Lfunc.init(Delta,hx,QUADRATIC_MODE);
      */

      /*
     if (IsZero(qi_pair<GF2EX>::get_rootD()*(qi_pair<GF2EX>::get_rootD() + hx) + Delta))
       return false;
      */
  }
  }


  //
  // QuadraticOrder<GF2EX>::IsEqual()
  //
  // Task:
  //      tests if the function fields are equal
  //

  template <>
  bool
  QuadraticOrder < GF2EX >::IsEqual (const QuadraticOrder < GF2EX > &QO) const
  {
    return ((Delta == QO.Delta) && (hx == QO.hx));
  }


  //
  // QuadraticOrder<GF2EX>::IsImaginary()
  //
  // Task:
  //      returns true if the function field is imaginary
  //
                                                                                
  template <> 
  bool 
  QuadraticOrder < GF2EX >::IsImaginary () const
  {
    return ((deg (Delta) & 1) && (deg(hx) <= ((deg(Delta)-1) >> 1)));
  }
                                                                                


  //
  // QuadraticOrder<GF2EX>::IsUnusual()
  //
  // Task:
  //      returns true if the function field is unusual
  //
                                                                                
  template <> 
  bool 
  QuadraticOrder < GF2EX >::IsUnusual () const
  {
    return (!(deg(Delta) & 1) && (deg(hx) == (deg(Delta) >> 1)) && !test_Dcoeff(Delta,hx));
  }
                                                                                


  //
  // QuadraticOrder<T>::IsReal()
  //
  // Task:
  //      returns true if the function field is real
  //

  template <>
  bool 
  QuadraticOrder < GF2EX >::IsReal () const
  {
    bool real1 = (!(deg(Delta) & 1) && (deg(hx) == (deg(Delta) >> 1)) && test_Dcoeff(Delta,hx));
    bool real2 = (deg(Delta) < 2*deg(hx));
    return (real1 || real2);
  }



  /*
  //
  // randomImaginaryOrder
  //
  // Task:
  //      generates a random imaginary order
  //

  template <>
    QuadraticOrder<GF2EX> & randomImaginaryOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    GF2EX disc;

    do
      {
	random (disc, 2 * size + 2);
	MakeMonic (disc);
      }
    while (deg (disc) != (2 * size + 1) || !DetIrredTest (disc));

    GF2EX newh;

    do
      {
	random (newh, size + 1);
	MakeMonic (newh);
      }
    while (deg (newh) < 1 || !DetIrredTest (newh));

    return new QuadraticOrder<GF2EX>(disc,newh);
  }



  //
  // randomUnusualOrder
  //
  // Task:
  //      generates a random real order
  //

  template <>
    QuadraticOrder<GF2EX> & randomUnusualOrderual (long size, bool prime)
  {
    // select a random discriminant with given genus

    GF2EX disc;
    GF2EX RD, newh;

    do
      {
	do
	  {
	    random (disc, 2 * size + 3);
	  }
	while (deg (disc) != (2 * size + 2) || !DetIrredTest (disc));

	do
	  {
	    random (newh, size + 2);
	    MakeMonic(newh);
	  }
	while (deg (newh) != size+1 || !DetIrredTest (hx));
      }
    while (!SqrRoot (RD, newh, disc) || test_Dcoeff(disc,newh));

    return new QuadraticOrder<GF2EX>(disc,newh);
  }



  //
  //
  // randomRealOrder
  //
  // Task:
  //      generates a random real order
  //

  template <>
    QuadraticOrder<GF2EX> & randomRealOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    GF2EX disc;
    GF2EX RD, newh;

    do
      {
	do
	  {
	    random (disc, 2 * size + 3);
	  }
	while (deg (disc) != (2 * size + 2) || !DetIrredTest (disc));

	do
	  {
	    random (newh, size + 2);
	    MakeMonic(newh);
	  }
	while (deg (newh) != size+1 || !DetIrredTest (hx));
      }
    while (!SqrRoot (RD, newh, disc));

    return new QuadraticOrder<GF2EX>(disc,newh);
  }
*/



/*
  //
  // QuadraticOrder<GF2EX>::use_Lfunction_tables(H)
  //
  // Task:
  //      initializes table-driven L-function approximations, assuming that
  //      discriminants will be less than H
  //
                                                                              
  template <> void QuadraticOrder < GF2EX >::use_Lfunction_tables (const GF2EX & H)
  {
                                                                              
  }
                                                                              


  //
  // generate_mu_table()
  //
  // Task:
  //      generates a table of mu values - the number of baby-steps per giant-step.
  //

  template <>
  void
  QuadraticOrder < GF2EX >::
  generate_mu_table (long size_start, long size_stop, long num)
  {
    register long i, j, size;

    qi_pair < GF2EX > A, B;
    timer tb, tg;
    RR rat, lD;

    for (i = size_start; i <= size_stop; ++i)
      {
	size = num / ((i / 2) + 1);

	randomize_real (i, true);

	qi_pair < GF2EX >::set_current_order (Delta, hx);


	A.assign_one ();
	tb.start_timer ();
	for (j = 0; j < num; ++j)
	  A.rho ();
	tb.stop_timer ();

	B = A;
	B.inverse_rho ();

	tg.start_timer ();
	for (j = 0; j < num; ++j)
	  nucomp (A, A, B);
	tg.stop_timer ();

	rat = to_RR (tg.user_time ()) / to_RR (tb.user_time ());
	lD = to_RR (g);

	cout << i << "   " << rat << flush;
	cout << "   " << rat / lD << flush;
	cout << "   " << get_mu (Delta) << flush;
	cout << endl;
      }
  }
  */



  //
  // operator >>
  //
  // Task:
  //      inputs a QuadraticOrder from the istream in.
  //

  template <>
  void
  QuadraticOrder<GF2EX>::read_from_file(std::istream & in)
  {
    GF2EX newfx, newhx;

    in >> newfx;
    in >> newhx;
    assign (newfx, newhx);
  }



  //
  // operator <<
  //
  // Task:
  //      outputs a QuadraticOrder to the ostream out.
  //

  template <>
  void
  QuadraticOrder<GF2EX>::write_to_file(std::ostream & out) const
  {
    out << "Quadratic order (genus " << g << "):" << endl;
    out << "      w = " << GF2E::modulus () << endl;
    out << "   f(x) = " << Delta << endl;
    out << "   h(x) = " << hx << endl;

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
  void QuadraticOrder <GF2EX>::bsgs_getentrysize (ZZ & entry_size, bool nodist)
  {
    entry_size = 3 * NTL_BITS_PER_LONG;	// space for pointers
    if (nodist)
      entry_size += 2 * g * NTL_BITS_PER_LONG * GF2E::WordLength ();
    else
      entry_size += 3 * g * NTL_BITS_PER_LONG * GF2E::WordLength ();
  }
  */

} // ANTL

