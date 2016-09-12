/**
 * @file QuadraticOrder_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from QuadraticOrder.hpp only.
 */

namespace ANTL
{
  //
  // constructor
  //     - sets resets all the class members which are vectors
  //

  template < class T > 
  QuadraticOrder < T>::QuadraticOrder (const T & D)
  {
    // test whether D is a valid discriminant
    if (deg (D) >! 3) {
      // assign new values
      Delta = D;
      clear (hx);

      // compute genus
      g = deg (Delta);
      if (g & 1)
        g = (g - 1) >> 1;
      else {
        g = (g >> 1) - 1;
	//        if (qi_pair<T>::get_rootD()*qi_pair<T>::get_rootD() == Delta)
	//return false;
      }

      /*
      // initialize invariant values
      clear (R);
      clear (h);
      CL.clear();
      gens.clear();
      clear (L);
  
      Lfunc.init(Delta,QUADRATIC_MODE);
      */
    }
  }



  //
  // destructor
  //

  template < class T > QuadraticOrder < T >::~QuadraticOrder ()
  {
  }



  /*
  //
  // QuadraticOrder<T>::verbose()
  //
  // Task:
  //      sets the verbosity of commands.  Currently, the following levels are
  //      supported:
  //         0 - nothing
  //         1 - run times and some run-time data (data only for subexp algs.)
  //         >1 - state information for subexponential algorithms
  //

  template < class T >
  void QuadraticOrder<T>::verbose (long state)
  {
    if (state <= 0)
      {
	info = 0;
      }
    else
      {
	info = state;
	L_function<T>::verbose(state);
      }
  }
  */



  //
  // QuadraticOrder<T>::IsEqual()
  //
  // Task:
  //      tests if the function fields are equal
  //

  template < class T >
  bool QuadraticOrder < T >::IsEqual (const QuadraticOrder < T > &QO) const
  {
    return (Delta == QO.Delta && hx == QO.hx);
  }



  //
  // operator ==
  //
  // Task:
  //      tests if QO1 and QO2 are equal (same discriminant)
  //

  template < class T > bool
  operator == (const QuadraticOrder < T > &QO1,
	       const QuadraticOrder < T > &QO2)
  {
    return (QO1.IsEqual (QO2));
  }



  //
  // operator !=
  //
  // Task:
  //      tests if QO1 and QO2 are not equal
  //

  template < class T > bool
  operator != (const QuadraticOrder < T > &QO1,
	       const QuadraticOrder < T > &QO2)
  {
    return (!QO1.IsEqual (QO2));
  }




  //
  // QuadraticOrder<T>::IsImaginary()
  //
  // Task:
  //      returns true if the function field is imaginary
  //

  template < class T > bool QuadraticOrder < T >::IsImaginary () const
  {
    return (deg (Delta) & 1);
  }



  //
  // QuadraticOrder<T>::IsUnusual()
  //
  // Task:
  //      returns true if the function field is real
  //

  template < class T > bool QuadraticOrder < T >::IsUnusual () const
  {
    return (!(deg(Delta) & 1) && !test_Dcoeff(Delta));
  }



  //
  // QuadraticOrder<T>::IsReal()
  //
  // Task:
  //      returns true if the function field is real
  //

  template < class T > bool QuadraticOrder < T >::IsReal () const
  {
    return (!(deg(Delta) & 1) && test_Dcoeff(Delta));
  }


  /*
  //
  // randomImaginaryOrder
  //
  // Task:
  //      generates a random imaginary order
  //

  template < class T >
  QuadraticOrder<T> & randomImaginaryOrder (long size, bool prime)
  {
    // select a random discriminant with given genus
    T disc;

    do
      {
	random (disc, 2 * size + 2);
	MakeMonic (disc);
      }
    while (deg (disc) != (2 * size + 1) || !DetIrredTest (disc));

    return new QuadraticOrder<T>(disc);
  }



  //
  // randomUnusualOrder
  //
  // Task:
  //      generates a random imaginary order
  //

  template < class T >
  QuadraticOrder<T> &
  randomUnusualOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    T disc;

    do
      {
	random (disc, 2 * size + 3);
	MakeMonic (disc);
      }
    while (deg (disc) != (2 * size + 2) || test_Dcoeff(disc));

    return new QuadraticOrder<T>(disc);
  }



  //
  // randomRealOrder
  //
  // Task:
  //      generates a random real order
  //

  template < class T >
  QuadraticOrder<T> &
  randomRealOrder (long size, bool prime)
  {
    // select a random discriminant with given genus

    T disc;

    do
      {
	random (disc, 2 * size + 3);
	MakeMonic (disc);
      }
    while (deg (disc) != (2 * size + 2) || !DetIrredTest (disc));

    return new QuadraticOrder<T>(disc);
  }
*/



  /*
  //
  // QuadraticOrder<T>::use_Lfunction_tables(H)
  //
  // Task:
  //      initializes table-driven L-function approximations, assuming that
  //      discriminants will be less than H
  //
                                                                              
  template < class T > 
  void QuadraticOrder < T >::use_Lfunction_tables (const T & H)
  {
  }
                                                                              
                                                                              
  //
  // QuadraticOrder<T>::set_Lfunction_global(H)
  //
  // Task:
  //      sets the global terms and precision constants for L_function
  //      computations
  //
                                                                              
  template < class T > 
  void QuadraticOrder < T >::set_Lfunction_global (const T & H)
  {
  }
                                                                              
      
  //
  // QuadraticOrder<T>::Lfunction
  //
  // Task:
  //      returns an approximation of L_K(1/q) computed via the analytic class
  //      number formula.
  //

  template < class T > RR QuadraticOrder < T >::Lfunction ()
  {
    if (!is_L_computed ())
      {
	// L has not been computed - compute it

	ZZ qg;

	power (qg, CARDINALITY<T>(), g);
	L = to_RR (class_number () * regulator ()) / to_RR (qg);
      }

    return L;
  }




  //
  // QuadraticOrder<T>::set_regulator(newR)
  //
  // Sets the regulator of the order - no error checking!
  //

  template < class T > void
  QuadraticOrder <
    T >::set_regulator (const ZZ & newR)
  {
    R = newR;
    Rbsgs = false;
    Rconditional = true;
  }



  //
  // QuadraticOrder<T>::set_class_number(newh)
  //
  // Sets the class number of the order - no error checking!
  //

  template < class T > void
  QuadraticOrder <
    T >::set_class_number (const ZZ & newh)
  {
    h = newh;
    hconditional = true;
  }

  template <class T> void QuadraticOrder<T>::use_bernstein(long bernnumb){
    rel_gen.use_bernstein(bernnumb);
  }



  //
  // QuadraticOrder<T>::exponent
  //
  // Task:
  //      computes the exponent of the class group.
  //

  template < class T > ZZ QuadraticOrder < T >::exponent ()
  {
    // compute class group, if needed
    if (!is_CL_computed ())
      class_group ();

    return CL.back();
  }



  //
  // QuadraticOrder<T>::p_rank
  //
  // Task:
  //      computes the p_rank of the class group.
  //

  template < class T >
  long QuadraticOrder<T>::p_rank (const ZZ & p)
  {
    long prank, i;
    ZZ temp;

    // compute class group, if needed
    if (!is_CL_computed ())
      class_group ();

    if (IsZero (h))
      prank = 0;
    else
      {
	// compute p-rank
	prank = 0;
	for (i=0; i<CL.size(); ++i) {
	  rem(temp,CL[i],p);
	  if (IsZero(temp))
	    ++prank;
	}
      }

    return prank;
  }






  // distance computation in BS-GS methods

  template < class T > void
  QuadraticOrder <
    T >::

  combine_BSGS (qo_distance < T > &dist, const qi_pair < T > &DD,
		const qo_hash_entry_real < T > *F)
  {
    dist = DD.get_distance () - F->get_d ();
  }


  template < class T > void
  QuadraticOrder <
    T >::

  combine_conj_BSGS (qo_distance < T > &dist, const qi_pair < T > &DD,
		     const qo_hash_entry_real < T > *F)
  {
    dist = DD.get_distance () + F->get_d () - deg (DD.get_a ());
  }


  template < class T > void
  QuadraticOrder <
    T >::

  combine_BSGS_easy (qo_distance < T > &dist, const qi_pair < T > &DD,
		     const qo_hash_entry_real < T > *F)
  {
    dist = F->get_d () - DD.get_distance ();
  }


  template < class T > void
  QuadraticOrder <
    T >::

  combine_conj_BSGS_easy (qo_distance < T > &dist, const qi_pair < T > &DD,
			  const qo_hash_entry_real < T > *F)
  {
    dist = R - F->get_d() + deg (DD.get_a ()) - DD.get_distance();
  }
  */



  /*
  //
  // read_from_file
  //
  // Task:
  //      inputs a QuadraticOrder from the istream in.
  //

  template < class T >
  void
  QuadraticOrder<T>::read_from_file(std::istream & in)
  {
    T newD;

    in >> newD;
    assign (newD);
  }
  */


  //
  // operator >>
  //
  // Task:
  //      inputs a QuadraticOrder from the istream in.
  //

  template < class T >
  std::istream & operator >> (std::istream & in, QuadraticOrder < T > &QO)
  {
    QO.read_from_file(in);
    return in;
  }


  //
  // operator <<
  //
  // Task:
  //      outputs a QuadraticOrder to the ostream out.
  //

  /*
  template < class T >
  void QuadraticOrder<T>::write_to_file(std::ostream & out) const
  {
    out << "Quadratic order (genus " << g << "):" << endl;
    output_FF (out, Delta);

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


  //
  //
  // operator <<
  //
  // Task:
  //      outputs a QuadraticOrder to the ostream out.
  //

  template < class T >
  std::ostream & operator << (std::ostream & out, const QuadraticOrder < T > &QO)
  {
    QO.write_to_file(out);
    return out;
  }

} // ANTL
