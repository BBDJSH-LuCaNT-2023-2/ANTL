/**
 * @file L_function.hpp
 * @author Leonard Nooy
 * @version $Header$
 *
 * Purpose: This file contains the class definition for the L_function
 *          class. The class includes routines to deal with specific
 *          cases of L functions, but can also be used to compute
 *          general cases of L functions.
 *          Note, this class outputs warnings to stdout, in the case
 *          a funtion is called and the mode, or template type is
 *          invalid.
 *
 * Revision History:
 *    May 16, 2003 - File created
 *    May 21, 2003 - Added more approximation functions, and some
 *                   character functions.
 *    May 23, 2003 - Discussed class with Mike, and refined some of
 *                   the function declarations. Added two new functions
 *                   approximateL1*
 *    May 27, 2003 - Added an instance of the Character class. This may
 *                   Change very quickly...
 *    May 28, 2003 - Added the function calculate_optimal_terms.
 *    June 3, 2003 - Added the get_poly_modq function for calculating the
 *                   Euler term of L(1) over function fields.
 *   June 24, 2003 - Added function to handle specific precision
 *                   calculations of the L1 over ZZ, and long
 *   July 22, 2003 - Added functions for calculating approximations to HR and h-
 *                   as well as cleaning up some code.
 *  August 6, 2003 - Added some set and get functions that Prof. Jacobson asked for
 *      Fall, 2003 - Replaced the float, double, xdouble, quad_float, and RR versions of
 *                   approximateL0_Quartic, approximateL0_ImaginaryNumberField,
 *                   approximateL0_RealNumberField, ApproximateL1, ApproximateL1_Quartic, and
 *                   approximateL1Bach_Quartic with templated friend functions (Rennie deGraaf)
 *      Fall, 2005 - modified to use improved, templated Character class, 
 *                   improved table-driven stuff, added long long
 *                   specialization (MJJ)
 */

#ifndef ANTL_L_FUNCTION_L_FUNCTION_H
#define ANTL_L_FUNCTION_L_FUNCTION_H

#include <cmath>
#include <NTL/ZZ.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEX.h>
#include <NTL/lzz_pEXFactoring.h>

#define L0_REF 0		/* array indexes into the nterms_used array */
#define L1_REF 1
#define L2_REF 2

#define QUADRATIC_L1_MODE 2
#define QUADRATIC_L0_MODE 3

#define MAX_L0_TABLE 40000

#include <ANTL/common.hpp>
#include <ANTL/Arithmetic/CC.hpp>
#include <ANTL/L_function/Character.hpp>
#include <ANTL/L_function/L_function_util.hpp>

namespace ANTL
{

  template < class T > class L_function
  {
  protected:
    T Delta;			/* Discriminant */
    T hx;				/* needed for GF2EX */
    int mode;			/* {quadratic, quartic} */
    long nterms_used[3];		/* contains the # of terms used in the last
					   calculation of L0, L1, L2. */

    CC < RR > L0_result;		/* stores previously calculated values for easy lookup. */
    RR L1_result;
    RR L2_result;

    int r2;                       /* number of linear factors of Delta (for FF cases) */

    static long info;

  public:
    Character<T> Chi;		/* this class contains the kronecker and quartic character stuff */

  protected:

    PrimeSeq primes;              /* list of small primes */
    RR PI;			/* global static value of PI set in the constructor */

    T a, b;			/* used to factor the Delta when in quartic L0 case, and Delta is a prime = 5 mod 8 */
    long epsilon;			/* stores the value of epsilon (1, -1) from the last quartic calculation */

    // global table and variables for table-driven L(1,X) calculations
    long table_size;              // bound on entries in tables
    long table_size_real;         // number of entries in tables

    double *log_table;            // table of i log i for Q <= i < 2Q
    double **L1table;             // log(P / (P - (D/P)), indexed by D mod P

    long *prime_list;             // list of small primes for L(1,X) approx
    long numP;                    // max prime in list

    // global table and variables for table-driven L(0,X) calculations
  public:
    short **chi_table;            // table of (D/n), indexed by D mod n

  protected:

    bool global_flag;		/* flag for usage of the global precision */
    long global_num_terms;
    long global_prec;

    void delete_tables();

    CC<RR> approximateL0_Quartic_default(double acc);
    RR approximateL0_ImaginaryNumberField_default (double acc);
    RR approximateL0_RealNumberField_default(double acc);
    CC<RR> approximateL0_default(long terms);
    CC<RR> approximateL0_default(double acc);

  public:
    /* CONSTRUCTOR / DESTRUCTOR */
    L_function ();
    ~L_function ();

    /* INITALIZATION FUNCTIONS */
    void init (const T & x, int m);
    void init (const T & x, const T & y, int m);
    void set_global_prec (const T & D, int m);
    void set_global_terms (const T & D, double acc, int m);

    static void verbose (long state);

    /* SET / GET FUNCTIONS */
    long get_epsilon (void) const;
    T get_a (void) const;		/* the two factors of the quartic discriminant */
    T get_b (void) const;
    int terms_used(int s) const {
      if (s > 2 || s < 0)  return 0;
      else return nterms_used[s];
    }

    /* GENERAL L FUNCTION */
    RR approximate (long s, long terms);	/* This function calls the more specific L1, L2, L function calls. */
    RR approximate (long s, double acc);

    RR approximateL (long s, long terms);

    /* L0 FUNCTIONS */
    long calculate_optimal_n_Quartic (const T & D, double acc);	/* used over ZZ, long, imaginary quartic discriminants */
    long calculate_precision_Quartic (const T & D);

    long calculate_optimal_n_RealNumberField (const T & D, double acc);	/* used over ZZ, long, Real Quadratic discriminants */
    long calculate_precision_RealNumberField (const T & D);

    long calculate_optimal_n_ImaginaryNumberField (const T & D, double acc);	/* used over ZZ, long, Imaginary Quadratic discriminants */
    long calculate_precision_ImaginaryNumberField (const T & D);

    CC < RR > approximateL0 (long terms);
    CC < RR > approximateL0 (double acc);

    CC < RR > approximateL0_Quartic (double acc);
    template <class TYPE> friend CC<RR> approximateL0_Quartic_impl(L_function<T>& lfunc, long terms);

    RR approximateL0_RealNumberField (double acc);
    template <class TYPE> friend RR approximateL0_RealNumberField_impl(L_function<T>& lfunc, long terms);

    RR approximateL0_ImaginaryNumberField (double acc);
    template <class TYPE> friend RR approximateL0_ImaginaryNumberField_impl(L_function<T>& lfunc, long terms);



    /* L1 FUNCTIONS */
    RR calculate_L1_error(const T & D, long Q);
    RR calculate_L1_error_ff(const T & D, long Q);
    long calculate_optimal_terms (const T & D, double acc);	/* used over ZZ, long, long long quadratic fields */
    long calculate_optimal_terms_new (const T & D, double acc);	/* used over ZZ, long, long long quadratic fields */
    long calculate_precision (const T & D);
    long calculate_optimal_terms_ff (double acc);    /* used over ZZ_pX, zz_pX, ZZ_pEX, zz_pEX, GF2EX */

    RR approximateL1 (long terms);
    RR approximateL1 (double acc);

    template <class TYPE> friend RR approximateL1_Quartic_impl(L_function<T>& lfunc, double acc);
    template <class TYPE> friend RR approximateL1Bach_Quartic_impl(L_function<T>& lfunc, long terms);

    template <class TYPE> friend RR approximateL1_impl(L_function<T>& lfunc, long terms);



    /* L2 FUNCTIONS */
    RR approximateL2 (long terms);

    /* MISC FUNCTIONS */
    RR calculate_phi_error (long terms);	/* used over zz_pX, ZZ_pX, etc.. */
    RR euler_term (long s, long n);

    // TABLE-DRIVEN FUNCTIONS
    void create_L0_tables (const T & D, double acc);
    void create_L1_tables (const T & D, double acc);

    RR approximateL0_RealNumberField_table (double acc);
    RR approximateL0_ImaginaryNumberField_table (double acc);

    RR approximateL1_table ();         // assumes global precision and terms set

  protected:
    template <class TYPE> friend RR approximateL0_RealNumberField_table_impl(L_function<T>& lfunc, long terms);
    template <class TYPE> friend RR approximateL0_ImaginaryNumberField_table_impl(L_function<T>& lfunc, long terms);
  };

  /* declarations of friend functions */
  /* ZZ specializations */
  template <class TYPE> CC <RR> approximateL0_Quartic_impl(L_function<ZZ>& lfunc, long terms);
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<ZZ>& lfunc, long terms);
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<ZZ>& lfunc, long terms);
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<ZZ>& lfunc, double acc);
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<ZZ>& lfunc, long terms);
  template <class TYPE> RR approximateL1_impl(L_function<ZZ>& lfunc, long terms);
  template <class TYPE> RR approximateL0_RealNumberField_table_impl(L_function<ZZ>& lfunc, long terms);
  template <class TYPE> RR approximateL0_ImaginaryNumberField_table_impl(L_function<ZZ>& lfunc, long terms);

  /* long specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<long>& lfunc, long terms);
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<long>& lfunc, double acc);
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<long>& lfunc, long terms);
  template <class TYPE> RR approximateL1_impl(L_function<long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_RealNumberField_table_impl(L_function<long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_ImaginaryNumberField_table_impl(L_function<long>& lfunc, long terms);

  /* long long specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<long long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<long long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<long long>& lfunc, long terms);
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<long long>& lfunc, double acc);
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<long long>& lfunc, long terms);
  template <class TYPE> RR approximateL1_impl(L_function<long long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_RealNumberField_table_impl(L_function<long long>& lfunc, long terms);
  template <class TYPE> RR approximateL0_ImaginaryNumberField_table_impl(L_function<long long>& lfunc, long terms);

  /* These specializations are not yet defined */
#define DEFAULT_APPROX_METHOD(TYPE, RET)				\
  cerr << "ApproximateL* method not defined for type " << #TYPE << endl; \
    return to<RET>(0.0);

  /* GF2EX specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<GF2EX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(GF2EX, CC<RR>) }
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<GF2EX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(GF2EX , RR) }
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<GF2EX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(GF2EX, RR) }
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<GF2EX>& lfunc, double acc) { DEFAULT_APPROX_METHOD(GF2EX, RR) }
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<GF2EX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(GF2EX, RR) }
  template <class TYPE> RR approximateL1_impl(L_function<GF2EX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(GF2EX, RR) }

  /* ZZ_pX specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<ZZ_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pX, CC<RR>) }
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<ZZ_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pX, RR) }
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<ZZ_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pX, RR) }
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<ZZ_pX>& lfunc, double acc) { DEFAULT_APPROX_METHOD(ZZ_pX, RR) }
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<ZZ_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pX, RR) }
  template <class TYPE> RR approximateL1_impl(L_function<ZZ_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pX, RR) }

  /* zz_pX specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<zz_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pX, CC<RR>) }
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<zz_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pX, RR) }
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<zz_pX>& lfunc, long terms) {DEFAULT_APPROX_METHOD(zz_pX, RR) }
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<zz_pX>& lfunc, double acc) { DEFAULT_APPROX_METHOD(zz_pX, RR) }
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<zz_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pX, RR) }
  template <class TYPE> RR approximateL1_impl(L_function<zz_pX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pX, RR) }

  /* ZZ_pEX specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<ZZ_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pEX, CC<RR>) }
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<ZZ_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pEX, RR) }
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<ZZ_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pEX, RR) }
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<ZZ_pEX>& lfunc, double acc) { DEFAULT_APPROX_METHOD(ZZ_pEX, RR) }
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<ZZ_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pEX, RR) }
  template <class TYPE> RR approximateL1_impl(L_function<ZZ_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(ZZ_pEX, RR) }

  /* zz_pEX specializations */
  template <class TYPE> CC<RR> approximateL0_Quartic_impl(L_function<zz_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pEX, CC<RR>) }
  template <class TYPE> RR approximateL0_RealNumberField_impl(L_function<zz_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pEX, RR) }
  template <class TYPE> RR approximateL0_ImaginaryNumberField_impl(L_function<zz_pEX>& lfunc, long terms) {DEFAULT_APPROX_METHOD(zz_pEX, RR) }
  template <class TYPE> RR approximateL1_Quartic_impl(L_function<zz_pEX>& lfunc, double acc) { DEFAULT_APPROX_METHOD(zz_pEX, RR) }
  template <class TYPE> RR approximateL1Bach_Quartic_impl(L_function<zz_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pEX, RR) }
  template <class TYPE> RR approximateL1_impl(L_function<zz_pEX>& lfunc, long terms) { DEFAULT_APPROX_METHOD(zz_pEX, RR) }

  //
  // Specialized Method Prototypes
  //
  // ??
  template <> RR L_function < GF2EX >::euler_term (long s, long n);
  template <> void L_function < ZZ >::init (const ZZ & x, int m);
  template <> CC < RR > L_function < ZZ >::approximateL0_Quartic (double acc);
  template <> RR L_function < ZZ >::approximateL0_ImaginaryNumberField (double acc);
  template <> RR L_function < ZZ >::approximateL0_RealNumberField (double acc);
  template <> CC < RR > L_function < ZZ >::approximateL0 (long terms);
  template <> CC < RR > L_function < ZZ >::approximateL0 (double acc);
  template <> RR L_function < ZZ >::approximateL1 (long n);
  template <> RR L_function < ZZ >::approximateL1 (double acc);
  template <> RR L_function < ZZ >::approximateL2 (long terms);
  template <> RR L_function < ZZ >::approximateL (long s, long terms);
  template <> RR L_function < ZZ >::approximate (long s, long n);
  template <> RR L_function < ZZ >::approximate (long s, double acc);
  //template <> void L_function < ZZ >::create_L1_tables (const ZZ & D, double acc);
  //template <> void L_function < ZZ >::create_L0_tables (const ZZ & D, double acc);
  template <> RR L_function < ZZ >::approximateL0_RealNumberField_table(double acc);
  template <> RR L_function < ZZ >::approximateL0_ImaginaryNumberField_table(double acc);
  template <> RR L_function < ZZ >::approximateL1_table ();

  template <> void L_function < long >::init (const long &x, int m);
  template <> CC < RR > L_function < long >::approximateL0_Quartic (double acc);
  template <> RR L_function < long >::approximateL0_ImaginaryNumberField (double acc);
  template <> RR L_function < long >::approximateL0_RealNumberField (double acc);
  template <> CC < RR > L_function < long >::approximateL0 (long terms);
  template <> CC < RR > L_function < long >::approximateL0 (double acc);
  template <> RR L_function < long >::approximateL1 (long n);
  template <> RR L_function < long >::approximateL1 (double acc);
  template <> RR L_function < long >::approximateL2 (long terms);
  template <> RR L_function < long >::approximateL (long s, long terms);
  template <> RR L_function < long >::approximate (long s, long n);
  template <> RR L_function < long >::approximate (long s, double acc);
  //template <> void L_function < long >::create_L1_tables (const long &D, double acc);
  //template <> void L_function < long >::create_L0_tables (const long &D, double acc);
  template <> RR L_function < long >::approximateL0_RealNumberField_table (double acc);
  template <> RR L_function < long >::approximateL0_ImaginaryNumberField_table(double acc);
  template <> RR L_function < long >::approximateL1_table ();
  template <> void L_function < long >::init (const long &x, int m);

  template <> CC < RR > L_function < long long >::approximateL0_Quartic (double acc);
  template <> RR L_function < long long >::approximateL0_ImaginaryNumberField (double acc);
  template <> RR L_function < long long >::approximateL0_RealNumberField (double acc);
  template <> CC < RR > L_function < long long >::approximateL0 (long terms);
  template <> CC < RR > L_function < long long >::approximateL0 (double acc);
  template <> RR L_function < long long >::approximateL1 (long n);
  template <> RR L_function < long long >::approximateL1 (double acc);
  template <> RR L_function < long long >::approximateL2 (long terms);
  template <> RR L_function < long long >::approximateL (long s, long terms);
  template <> RR L_function < long long >::approximate (long s, long n);
  template <> RR L_function < long long >::approximate (long s, double acc);
  template <> RR L_function < long long >::approximateL0_RealNumberField_table (double acc);
  template <> RR L_function < long long >::approximateL0_ImaginaryNumberField_table(double acc);
  template <> RR L_function < long long >::approximateL1_table ();
  template <> void L_function < long long >::init (const long long &x, int m);

} // ANTL

/* base instantiation, ZZ_pX, zz_pX, etc... */
#include "../../../src/L_function/L_function_impl.hpp" 
#include "../../../src/L_function/L_function_FF_impl.hpp"

#endif // guard
