/**
 * @file QuadraticOrder_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from QuadraticOrder.hpp only.
 */

#include <iostream>

namespace ANTL
{
  //
  // constructor
  //     - sets resets all the class members which are vectors
  //

  template < class T > QuadraticOrder<T>::QuadraticOrder (const T & D)
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
      */
      Lfunc.init(Delta,QUADRATIC_MODE);
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

// Multiply Strategy Setters
  template <class T>
  void QuadraticOrder<T>::set_red_best(ReduceStrategy<T> &A) {
    A.init(Delta, hx, g);
    red_best = &A;
  }

  template <class T>
  void QuadraticOrder<T>::set_red_plain_imag(ReducePlainImag<T> &A) {
    A.init(Delta, hx, g);
    red_plain_imag = &A;

    //If red_best is not set, or if Plain is the preferred strategy, set red_best to A as well
    if(red_best == nullptr || (preferred_red == 0 && this->IsImaginary())) {
      red_best = &A;
    }
  }

  template <class T>
  void QuadraticOrder<T>::set_red_plain_real(ReducePlainReal<T> &A) {
    A.init(Delta, hx, g);
    red_plain_real = &A;

    //If red_best is not set, or if Plain is the preferred strategy, set red_best to A as well
    if(red_best == nullptr || (preferred_red == 0 && this->IsReal())) {
      red_best = &A;
    }
  }

  template <class T>
  void QuadraticOrder<T>::set_red_fast(ReduceFast<T> &A) {
    A.init(Delta, hx, g);
    red_fast = &A;

    //If red_best is not set, or if Fast is the preferred strategy, set red_best to A as well
    if(red_best == nullptr || preferred_red == 1) {
      red_best = &A;
    }
  }

// Multiply Strategy Setters
  template <class T> void QuadraticOrder<T>::set_mul_best(MultiplyStrategy<T> &A) {
    mul_best = &A;
  }

  template <class T> void QuadraticOrder<T>::set_mul_plain(MultiplyPlain<T> &A) {
    mul_plain = &A;

    //If mul_best is not set, or if Plain is the preferred strategy, set mul_best to A as well
    if(mul_best == nullptr || preferred_mul == 0) {
      mul_best = &A;
    }
  }

  template <class T> void QuadraticOrder<T>::set_mul_nucomp(MultiplyNucomp<T> &A) {
    mul_nucomp = &A;

    //If mul_best is not set, or if Nucomp is the preferred strategy, set mul_best to A as well
    if(mul_best == nullptr || preferred_mul == 1) {
      mul_best = &A;
    }
  }

// Square Strategy Setters
  template <class T> void QuadraticOrder<T>::set_sqr_best(SquareStrategy<T> &A) {
    sqr_best = &A;
  }

  template <class T> void QuadraticOrder<T>::set_sqr_plain(SquarePlain<T> &A) {
    sqr_plain = &A;

    //If sqr_best is not set, or if Plain is the preferred strategy, set sqr_best to A as well
    if(sqr_best == nullptr || preferred_sqr == 0) {
      sqr_best = &A;
    }
  }

  template <class T> void QuadraticOrder<T>::set_sqr_nudupl(SquareNudupl<T> &A) {
    sqr_nudupl = &A;

    //If sqr_best is not set, or if Nudupl is the preferred strategy, set sqr_best to A as well
    if(sqr_best == nullptr || preferred_sqr == 0) {
      sqr_best = &A;
    }
  }

// Cube Strategy Setters
  template <class T> void QuadraticOrder<T>::set_cube_best(CubeStrategy<T> &A) {
    cube_best = &A;
  }

  template <class T> void QuadraticOrder<T>::set_cube_plain(CubePlain<T> &A) {
    cube_plain = &A;

    //If cube_best is not set, or if Plain is the preferred strategy, set cube_best to A as well
    if(cube_best == nullptr || preferred_cube == 0) {
      cube_best = &A;
    }
  }

  template <class T> void QuadraticOrder<T>::set_cube_nucube(CubeNucube<T> &A) {
    cube_nucube = &A;

    //If cube_best is not set, or if Nucube is the preferred strategy, set cube_best to A as well
    if(cube_best == nullptr || preferred_cube == 1) {
      cube_best = &A;
    }
  }

  template <class T> void QuadraticOrder<T>::set_cube_mulsqr(CubeMulSqr<T> &A) {
    cube_mulsqr = &A;

    //If cube_best is not set, or if Mulsqr is the preferred strategy, set cube_best to A as well
    if(cube_best == nullptr || preferred_cube == 2) {
      cube_best = &A;
    }
  }

  template <class T> ReduceStrategy<T> *   QuadraticOrder<T>::get_red_best()       {return red_best;}
  template <class T> ReducePlainImag<T> *  QuadraticOrder<T>::get_red_plain_imag() {return red_plain_imag;}
  template <class T> ReducePlainReal<T> *  QuadraticOrder<T>::get_red_plain_real() {return red_plain_real;}
  template <class T> ReduceFast<T> *       QuadraticOrder<T>::get_red_fast()       {return red_fast;}

  template <class T> MultiplyStrategy<T> * QuadraticOrder<T>::get_mul_best()      {return mul_best;}
  template <class T> MultiplyPlain<T> *    QuadraticOrder<T>::get_mul_plain()     {return mul_plain;}
  template <class T> MultiplyNucomp<T> *   QuadraticOrder<T>::get_mul_nucomp()    {return mul_nucomp;}

  template <class T> SquareStrategy<T> *   QuadraticOrder<T>::get_sqr_best()      {return sqr_best;}
  template <class T> SquarePlain<T> *      QuadraticOrder<T>::get_sqr_plain()     {return sqr_plain;}
  template <class T> SquareNudupl<T> *     QuadraticOrder<T>::get_sqr_nudupl()    {return sqr_nudupl;}

  template <class T> CubeStrategy<T> *     QuadraticOrder<T>::get_cube_best()     {return cube_best;}
  template <class T> CubePlain<T> *        QuadraticOrder<T>::get_cube_plain()    {return cube_plain;}
  template <class T> CubeNucube<T> *       QuadraticOrder<T>::get_cube_nucube()   {return cube_nucube;}
  template <class T> CubeMulSqr<T> *       QuadraticOrder<T>::get_cube_mulsqr()   {return cube_mulsqr;}


  template < class T > long QuadraticOrder<T>::bsgs_getl(const ZZ & K, ZZ & N, ZZ & entry_size, RR & mu, bool nodist) {
    ZZ maxN, rootK;
    RR n, temp;
    long
    l;

    // get mu (# baby steps per giant step)
    mu = get_mu (Delta);

    // n = number of slave processes
//     if (parallel)
//       n = parallel;
//     else
      set (n);

    // compute max number of baby steps to store
    bsgs_getentrysize (entry_size, nodist);
    maxN = 1 + (MAXMEMORY / (to_ZZ (n) * entry_size));

    // compute number of baby steps assuming l=1, N=sqrt(KG/2n)
    temp = floor (SqrRoot (to_RR (K) * mu / (to_RR (1) * n)));
    conv (N, temp);
    l = 1;

    if (N > maxN) {
      // compute l = sqrt(KG/2n) / N;
      //    N = maxN;
      entry_size *= 3;
      entry_size >>= 1;
      N = 1 + (MAXMEMORY / (to_ZZ (n) * entry_size));
      temp = floor ((SqrRoot (to_RR (K) * mu / (to_RR (2) * n)) / to_RR (N)));
      conv (l, temp);
      if (l < 1) {
        l = 1;
      }
    }

    return l;
  }
// QuadraticOrder<T>::approximate_hR()
// Task: returns an approximation of hR

  template < class T > ZZ QuadraticOrder< T >::approximate_hR() {
    ZZ hR, temp, q = CARDINALITY<T>();
    long n;

    if (g == 1) {
      hR = q;
      ++hR;
      return hR;
    }

    if (g == 2) {
      hR = q * q + 6 * q + 1;
      return hR;
    }

    n = 2*g - 1;

    if ((g % 5) == 2) {
      n /= 5;
    }
    else {
      double tn = (double) n;
      n = (long) ::round(tn/5.0);
    }

    if (q < 5) {
      if (n < 8) {
        n = 8;
      }
    }
    else if (q < 7) {
      if (n < 4) {
        n = 4;
      }
    }
    else if (q < 11) {
      if (n < 3) {
        n = 3;
      }
    }
    else if (q < 13) {
      if (n < 2) {
        n = 2;
      }
    }

    RR FI = Lfunc.approximateL1 (n);

    hR = CeilToZZ(FI);

    return hR;
  }

// QuadraticOrder<T>::lower_bound_hR()
// Task: returns a lower bound of hR such that L < hR < 2L

  template < class T > ZZ QuadraticOrder< T >::lower_bound_hR() {
    ZZ hR, q = CARDINALITY<T>();
    RR temp, temp2, rootq = SqrRoot (to_RR (q));

    // compute temp = (sqrt(q) + 1)^2g - (sqrt(q) - 1)^2g
    temp2 = rootq;
    temp = temp2;
    ++temp;
    --temp2;
    power(temp,temp,g << 1);
    power(temp2,temp2,g << 1);
    temp -= temp2;

    // if width of Hasse interval is <= 2, use lower Hasse bound
    if (temp <= to_RR(2)) {
      hR = CeilToZZ (temp2);

      return hR;
    }

    // approximate L1 to accuracy log(sqrt(2))
    double acc = log(sqrt((double) 2));
    RR FI = Lfunc.approximateL1 (acc);

    // lower bound is FI/sqrt(2)
    hR = CeilToZZ (FI / SqrRoot (to_RR (2)));

    return hR;
  }

// QuadraticOrder<T>::estimate_hR_error(long n)
// Task: returns L such that |hR - hR'| < exp(L)^2

  template < class T > ZZ QuadraticOrder< T >::estimate_hR_error() {
    ZZ err, qn, q = CARDINALITY<T>();
    long n = Lfunc.terms_used(1);

    if ((g == 1) && (n == 0)) {
      SqrRoot (err, q);
      ++err;
      err <<= 1;
      return err;
    }

    if ((g == 2) && (n == 0)) {
      SqrRoot (err, q);
      ++err;
      err *= ((q + 1) << 2);
      return err;
    }

    err = CeilToZZ(Lfunc.calculate_L1_error_ff(Delta,n));

    return err;
  }

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
