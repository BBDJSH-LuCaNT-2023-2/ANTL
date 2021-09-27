/**
 * @file quadratic_form.hpp
 * @author Michael Jacobson
 * @version $Header$
 */

#ifndef ANTL_QUADRATIC_QUADRATIC_FORM_H
#define ANTL_QUADRATIC_QUADRATIC_FORM_H

#include <ANTL/common.hpp>

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>

namespace ANTL
{
  // Forward declarations
  template <class T> class quadratic_form;
  template <class T> class quadratic_number;
  template <class T> class qo_hash_entry;
  template <class T> class qo_hash_entry_small;
  template <class T, class S> class qo_hash_entry_int;
  template <class T> class qo_hash_entry_vec;

  //
  // Declare Friends
  //

  template < class T >
  void nugcd(T & x, T & b, T & y, T & a);

  template < class T >
  void nugcd(T & x, T & b, T & y, T & a, bool & flag);

  template < class T >
  void nugcd_lehmer(T & x, T & b, T & y, T & a, bool & flag);

  template < class T >
  void nugcd_lehmer2(T & x, T & b, T & y, T & a, bool & flag);

  template < class T >
  void nugcd_lehmer3(T & x, T & b, T & y, T & a, T & saved, bool & flag);

  template < class T >
  T multiply_base(quadratic_form<T> & C,
		  const quadratic_form<T> & A,
		  const quadratic_form<T> & B);

  template < class T >
  T square_base(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void nucomp_base(quadratic_form<T> & C, const quadratic_form<T> & A,
		   const quadratic_form<T> & B);

  template < class T >
  void enucomp_base(quadratic_form<T> & C, T & x, T & y, T & z,
		    const quadratic_form<T> & A, const quadratic_form<T> & B);

  template < class T >
  void nudupl_base (quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  T nucomp_base (quadratic_form < T > &C, const quadratic_form < T > &A, const quadratic_form < T > &B, T & G, T & D);

  template < class T >
  T nudupl_base (quadratic_form < T > &C, const quadratic_form < T > &A, T & G, T & B);

  template < class T >
  quadratic_form<T> conjugate(const quadratic_form<T> & A);

  template < class T >
  void multiply_imag(quadratic_form<T> & C, const quadratic_form<T> & A,
		     const quadratic_form<T> & B);

  template < class T >
  void multiply_real(quadratic_form<T> & C, const quadratic_form<T> & A,
		     const quadratic_form<T> & B);

  template < class T >
  void multiply_real(quadratic_form<T> & C, const quadratic_form<T> & A,
		     const quadratic_form<T> & B, quadratic_number<T> & gamma,
		     const quadratic_number<T> &alpha,
		     const quadratic_number<T> &beta);


  template < class T >
  void multiply(quadratic_form<T> & C, const quadratic_form<T> & A,
		const quadratic_form<T> & B);

  template < class T >
  void square_imag(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void square_real(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void square(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void power(quadratic_form<T> & C, const quadratic_form<T> & A, const ZZ & n);

  // WHY IS THIS HERE?   MJJ
  template < class T >
  void power_reduce(quadratic_form<T> & C, const quadratic_form<T> & A, const ZZ & n);

  template < class T >
  void nucomp_imag(quadratic_form<T> & C, const quadratic_form<T> & A,
		   const quadratic_form<T> & B);

  template < class T >
  void nucomp_real(quadratic_form<T> & C, const quadratic_form<T> & A,
		   const quadratic_form<T> & B);

  template < class T >
  void enucomp_real(quadratic_form<T> & C, T & x, T & y, T & z,
		    const quadratic_form<T> & A, const quadratic_form<T> & B);

  template < class T >
  void nucomp(quadratic_form<T> & C, const quadratic_form<T> & A,
	      const quadratic_form<T> & B);

  template < class T >
  void nudupl_imag(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void nudupl_real(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void nudupl(quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void nupower(quadratic_form<T> & C, const quadratic_form<T> & A, const ZZ & n);

  template < class T >
  quadratic_form<T> operator - (const quadratic_form<T> & A);

  template < class T >
  quadratic_form<T> operator * (const quadratic_form<T> & A,
				const quadratic_form<T> & B);

  template < class T >
  void swap(quadratic_form<T> & A, quadratic_form<T> & B);

  template < class T >
  bool operator == (const quadratic_form<T> & A, const quadratic_form<T> & B);

  template < class T >
  bool operator != (const quadratic_form<T> & A, const quadratic_form<T> & B);

  template < class T >
  std::istream & operator >> (std::istream & in, quadratic_form<T> & A);

  template < class T >
  std::ostream & operator << (std::ostream & out, const quadratic_form<T> & A);

  template < class T >
  void nucompcd_base(quadratic_form<T> & C, const quadratic_form<T> & A,
		     const quadratic_form<T> & B);

  template < class T >
  void nuduplcd_base (quadratic_form<T> & C, const quadratic_form<T> & A);

  template < class T >
  void nupower_cd(quadratic_form<T> & C, const quadratic_form<T> & A, const ZZ & n);

  //
  // Class: quadratic_form<T>
  //
  // This class represents an element in the class group of a real quadratic
  //    order, i.e., a reduced, primitive, invertible ideal
  //    aZ + (b+sqrt(Delta))/2 Z  where Delta is the discriminant of the real
  //    quadratic order.  This class is derived from qi_class and shares all
  //    its properties.  In addition, a distance is associated with each
  //    instance of this class, so that computations in the infrastructure are
  //    possible.
  //

  template < class  T >
  class quadratic_form
  {
  protected:

    // ideal coefficients
    T a;
    T b;
    T c;


    // current quadratic_order: discriminant, radicand, and exponent bound
    static T hx;  // only for GF2EX
    static T Delta;
    static T rootD;

    static T D4;
    static long D4_long;
    static long genus;
    static short order_type;

  public:
#ifdef DEGS
    static ZZ total_comp;
    static ZZ total_nc;

    static ZZ max_red;
    static ZZ total_red;
    static ZZ avg_red;
    static ZZ max_red_nc;
    static ZZ total_red_nc;
    static ZZ avg_red_nc;

    static ZZ nucomp_G1;
    static ZZ nudupl_G1;
#endif

  protected:
    void test_form(const char *);

    // reduction routines
    void normalize_imag_mult();
    void normalize_imag();
    void normalize_real_mult();
    void normalize_real();
    void normalize_mult()
    {
      if (is_real()) normalize_real_mult ();
      else normalize_imag_mult ();
    }

    void normalize_unusual ();
    void normalize_unusual (quadratic_number<T> & alpha);


    friend void nugcd <T> (T & x, T & b, T & y, T & a);
    friend void nugcd <T> (T & x, T & b, T & y, T & a, bool & flag);
    friend void nugcd_lehmer <T> (T & x, T & b, T & y, T & a, bool & flag);
    friend void nugcd_lehmer2 <T> (T & x, T & b, T & y, T & a, bool & flag);
    friend void nugcd_lehmer3 <T> (T & x, T & b, T & y, T & a, T & saved, bool & flag);




    friend T multiply_base <> (quadratic_form<T> & C,
			       const quadratic_form<T> & A,
			       const quadratic_form<T> & B);
    friend T square_base <T> (quadratic_form<T> & C,
			      const quadratic_form<T> & A);

    friend void nucomp_base <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A,
				 const quadratic_form<T> & B);

    friend void enucomp_base <T> (quadratic_form<T> & C,
				  T & x, T & y, T & z,
				  const quadratic_form<T> & A,
				  const quadratic_form<T> & B);

    friend void nudupl_base <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A);

    friend T nucomp_base < T > (quadratic_form < T > &C,
				const quadratic_form < T > &A,
				const quadratic_form < T > &B, T & G, T & BB);

    friend T nudupl_base < T > (quadratic_form < T > &C,
				const quadratic_form < T > &A, T & G, T & B);

    friend void nucompcd_base <T> (quadratic_form<T> & C,
				   const quadratic_form<T> & A,
				   const quadratic_form<T> & B);

    friend void nuduplcd_base <T> (quadratic_form<T> & C,
				   const quadratic_form<T> & A);

  public:

    //
    // constructors and destructor
    //

    quadratic_form();
    quadratic_form(const quadratic_form<T> & A);
    quadratic_form(const qo_hash_entry<T> & A);


    ~quadratic_form();



    //
    // initialization
    //

    static void set_current_order(const T & newf, const T & newh);
    static void set_current_order(const T & newDelta);


    //
    // assignment
    //

    void assign_one();
    bool assign_prime(const T & p);
    bool assign(const T & a2, const T & b2);

    //bool assign(const long & a2, const long & b2) // FIXME-JLM
    //{ return assign(to_ZZ(a2), to_ZZ(b2)); }
    // WHY IS THIS HERE?  MJJ
    bool assign(const long& a2, const ZZ & b2)
    {return assign(to_ZZ(a2), b2); }

    void assign (const qo_hash_entry<T> & B);
    void assign(const quadratic_form<T> & B);
    quadratic_form<T> & operator = (const quadratic_form<T> & A);



    //
    // access functions
    //

    static T get_hx()  { return hx; };
    static T sqrtD()  { return rootD; };
    static T discriminant()  { return Delta; };
    static long get_genus() { return genus; };
    static short order_signature();
    static bool is_real();
    static bool is_imaginary();
    static bool is_unusual();

    T get_a() const  { return a; };
    T get_b() const  { return b; };
    T get_c() const  { return c; };
    T content () const  { return GCD(a,GCD(b,c)); };


    //
    // arithmetic operations
    //

    void normalize ();

    void reduce_imag_mult();
    void reduce_imag_nc ();
    void reduce_imag();
    void reduce_real();
    void reduce();

    void reduce_imag(quadratic_number<T> & alpha);
    void reduce_real(quadratic_number<T> & alpha);
    void reduce_real(vector< quadratic_number<T> > & alpha);
    void reduce(quadratic_number<T> & alpha);

    void rho();
    void inverse_rho();
    void rho(quadratic_number<T> & alpha);
    void inverse_rho(quadratic_number<T> & alpha);

    friend quadratic_form<T> conjugate <T> (const quadratic_form<T> & A);

    friend void multiply_imag <T> (quadratic_form<T> & C,
				   const quadratic_form<T> & A,
				   const quadratic_form<T> & B);

    friend void multiply_real <T> (quadratic_form<T> & C,
				   const quadratic_form<T> & A,
				   const quadratic_form<T> & B);

    friend void multiply_real <T> (quadratic_form<T> & C,
				   const quadratic_form<T> & A,
				   const quadratic_form<T> & B,
				   quadratic_number<T> & gamma,
				   const quadratic_number<T> & alpha,
				   const quadratic_number<T> & beta);


    friend void multiply <T> (quadratic_form<T> & C,
			      const quadratic_form<T> & A,
			      const quadratic_form<T> & B);

    friend void square_imag <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A);
    friend void square_real <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A);
    friend void square <T> (quadratic_form<T> & C,
			    const quadratic_form<T> & A);

    friend void power <T> (quadratic_form<T> & C,
			   const quadratic_form<T> & A,
			   const ZZ & n);

    // WHY IS THIS HERE?  MJJ
    friend void power_reduce <T> (quadratic_form<T> & C,
				  const quadratic_form<T> & A,
				  const ZZ & n);

    friend void nucomp_imag <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A,
				 const quadratic_form<T> & B);
    friend void nucomp_real <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A,
				 const quadratic_form<T> & B);
    friend void enucomp_real <T> (quadratic_form<T> & C,
				  T & x, T & y, T & z,
				  const quadratic_form<T> & A,
				  const quadratic_form<T> & B);
    friend void nucomp <T> (quadratic_form<T> & C,
			    const quadratic_form<T> & A,
			    const quadratic_form<T> & B);

    friend void nudupl_imag <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A);
    friend void nudupl_real <T> (quadratic_form<T> & C,
				 const quadratic_form<T> & A);
    friend void nudupl <T> (quadratic_form<T> & C,
			    const quadratic_form<T> & A);

    friend void nupower <T> (quadratic_form<T> & C,
			     const quadratic_form<T> & A,
			     const ZZ & n);

    friend void nupower_cd <T> (quadratic_form<T> & C,
				const quadratic_form<T> & A,
				const ZZ & n);

    friend quadratic_form<T> operator - <T>  (const quadratic_form<T> & A);
    friend quadratic_form<T> operator * <T> (const quadratic_form<T> & A,
					     const quadratic_form<T> & B);
    quadratic_form<T> & operator *= (const quadratic_form<T> & A);



    //
    // basic functions
    //

    friend void swap <T> (quadratic_form<T> & A, quadratic_form<T> & B);
    qo_hash_entry<T> hash() const;
    qo_hash_entry_small<T> hash_small() const;

    template < class S >
    qo_hash_entry_int<T,S> hash_int(const S & newd) const;

    qo_hash_entry_vec<T> hash_vec(const vec_ZZ & newd) const;



    //
    // comparisons
    //

    bool is_one() const;
    bool is_equal(const quadratic_form<T> & B) const;
    bool is_equal(const qo_hash_entry<T> & B) const;
    bool is_equal(const qo_hash_entry_small<T> & B) const;

    friend bool operator == <T> (const quadratic_form<T> & A,
				 const quadratic_form<T> & B);
    friend bool operator != <T> (const quadratic_form<T> & A,
				 const quadratic_form<T> & B);


    //
    // input/output
    //

    friend std::istream & operator >> <T> (std::istream & in, quadratic_form<T> & A);
    friend std::ostream & operator << <T> (std::ostream & out, const quadratic_form<T> & A);
  };

  //
  // Declare Specialized Methods
  //

  template <> void quadratic_form<GF2EX>::test_form(const char * msg);
  template <> void quadratic_form<GF2EX>::normalize_imag();
  template <> void quadratic_form<GF2EX>::normalize_real_mult();
  template <> void quadratic_form<GF2EX>::normalize_real();
  template <> void quadratic_form<GF2EX>::normalize_unusual();
  template <> void quadratic_form<GF2EX>::normalize_unusual(quadratic_number<GF2EX> & alpha);
  template <> void quadratic_form<GF2EX>::reduce_imag_mult();
  template <> void quadratic_form<GF2EX>::reduce_imag_nc();
  template <> void quadratic_form<GF2EX>::reduce_real();
  template <> void quadratic_form<GF2EX>::reduce_real(quadratic_number<GF2EX> & alpha);
  template <> void quadratic_form<GF2EX>::rho();
  template <> void quadratic_form<GF2EX>::inverse_rho();
  template <> void quadratic_form<GF2EX>::rho(quadratic_number<GF2EX> & alpha);
  template <> void quadratic_form<GF2EX>::inverse_rho(quadratic_number<GF2EX> & alpha);
  template <> GF2EX multiply_base(quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B);
  template <> GF2EX square_base<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A);
  template <> void nugcd<GF2EX> (GF2EX & x, GF2EX & b, GF2EX & y, GF2EX & a);
  template <> void nucomp_base <GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B);
  template <> void nudupl_base<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A);
  template <> GF2EX nucomp_base <GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B, GF2EX & G, GF2EX & BB);
  template <> GF2EX nudupl_base<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, GF2EX & G, GF2EX & B);
  template <> quadratic_form<GF2EX>::quadratic_form(const qo_hash_entry<GF2EX> & A);
  template <> void quadratic_form<GF2EX>::set_current_order(const GF2EX & newf, const GF2EX & newh);
  template <> short quadratic_form<GF2EX>::order_signature();
  template <> void quadratic_form<GF2EX>::assign_one();
  template <> void quadratic_form <GF2EX>::assign (const qo_hash_entry < GF2EX > &B);
  template <> bool quadratic_form<GF2EX>::assign_prime (const GF2EX & p);
  template <> bool quadratic_form<GF2EX>::assign(const GF2EX & a2, const GF2EX & b2);
  template <> quadratic_form<GF2EX> conjugate<GF2EX> (const quadratic_form<GF2EX> & A);
  template <> void nucompcd_base <GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B);
  template <> void nuduplcd_base<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A);
  /*
    template <> void nucomp_imag<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B);
    template <> void nucomp_real<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B);
    template <> void nucomp<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A, const quadratic_form<GF2EX> & B);
    template <> void nudupl_imag<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A);
    template <> inline void nudupl_real<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A);
    template <> void nudupl<GF2EX> (quadratic_form<GF2EX> & C, const quadratic_form<GF2EX> & A);
  */

  template <> void quadratic_form<ZZ>::test_form(const char * msg);
  template <> void quadratic_form<ZZ>::normalize_imag_mult();
  template <> void quadratic_form<ZZ>::normalize_real_mult();
  template <> void quadratic_form<ZZ>::normalize_imag();
  template <> void quadratic_form<ZZ>::normalize_real();
  template <> void quadratic_form<ZZ>::normalize_unusual();
  template <> void quadratic_form<ZZ>::normalize_unusual(quadratic_number<ZZ> & alpha);
  template <> void quadratic_form<ZZ>::reduce_imag_mult();
  template <> void quadratic_form<ZZ>::reduce_imag_nc();
  template <> void quadratic_form<ZZ>::reduce_real();
  template <> void quadratic_form<ZZ>::reduce_real(quadratic_number<ZZ> &alpha);
  template <> void quadratic_form<ZZ>::reduce_real(vector< quadratic_number<ZZ> > &alpha);
  template <> void quadratic_form<ZZ>::rho();
  template <> void quadratic_form<ZZ>::inverse_rho();
  template <> void quadratic_form<ZZ>::rho(quadratic_number<ZZ> & alpha);
  template <> void quadratic_form<ZZ>::inverse_rho(quadratic_number<ZZ> & alpha);
  template <> ZZ multiply_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A, const quadratic_form<ZZ> & B);
  template <> ZZ square_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A);
  // REMOVE?  MJJ	template <> void nugcd<ZZ> (ZZ & x, ZZ & b, ZZ & y, ZZ & a);
  template <> void nugcd<ZZ> (ZZ & x, ZZ & b, ZZ & y, ZZ & a, bool & flag);
  template <> void nugcd_lehmer<ZZ> (ZZ & x, ZZ & v, ZZ & y, ZZ & u, bool & flag);
  template <> void nugcd_lehmer2<ZZ> (ZZ & x, ZZ & v, ZZ & y, ZZ & u, bool & flag);
  template <> void nugcd_lehmer3<ZZ> (ZZ & x, ZZ & b, ZZ & y, ZZ & a, ZZ & saved, bool & flag);

  template <> void nucomp_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A, const quadratic_form<ZZ> & B);
  template <> void enucomp_base<ZZ> (quadratic_form<ZZ> & C, ZZ & x, ZZ & y, ZZ & z, const quadratic_form<ZZ> & A, const quadratic_form<ZZ> & B);
  template <> void nudupl_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A);
  template <> ZZ nucomp_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A, const quadratic_form<ZZ> & B, ZZ & G, ZZ & BB);
  template <> ZZ nudupl_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A, ZZ & G, ZZ & B);
  template <> void quadratic_form<ZZ>::set_current_order(const ZZ & newDelta);
  template <> void quadratic_form<ZZ>::set_current_order(const ZZ & newDelta, const ZZ & newh);
  template <> short quadratic_form<ZZ>::order_signature();
  template <> quadratic_form<ZZ>::quadratic_form(const qo_hash_entry<ZZ> & A);
  template <> void quadratic_form<ZZ>::assign_one();
  template <> bool quadratic_form<ZZ>::assign_prime (const ZZ & p);
  template <> bool quadratic_form<ZZ>::assign(const ZZ & a2, const ZZ & b2);
  template <> bool quadratic_form<ZZ>::is_equal(const quadratic_form<ZZ> & B) const;
  template <> bool quadratic_form<ZZ>::is_equal(const qo_hash_entry<ZZ> & B) const;
  template <> void nucompcd_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A, const quadratic_form<ZZ> & B);
  template <> void nuduplcd_base<ZZ> (quadratic_form<ZZ> & C, const quadratic_form<ZZ> & A);

  template <> void quadratic_form < long >::test_form (const char *msg);
  template <> void quadratic_form < long >::normalize_imag_mult ();
  template <> void quadratic_form < long >::normalize_imag ();
  template <> void quadratic_form < long >::normalize_real_mult ();
  template <> void quadratic_form < long >::normalize_real ();
  template <> void quadratic_form<long>::normalize_unusual();
  template <> void quadratic_form<long>::normalize_unusual(quadratic_number<long> & alpha);
  template <> void quadratic_form < long >::rho ();
  template <> void quadratic_form < long >::inverse_rho ();
  template <> void quadratic_form < long >::rho (quadratic_number<long> & alpha);
  template <> void quadratic_form < long >::inverse_rho (quadratic_number<long> & alpha);
  template <> void quadratic_form < long >::reduce_imag_mult ();
  template <> void quadratic_form < long >::reduce_imag_nc ();
  template <> void quadratic_form < long >::reduce_real ();
  template <> long multiply_base < long >(quadratic_form < long >&C, const quadratic_form < long >&A, const quadratic_form < long >&B);
  template <> long square_base < long >(quadratic_form < long >&C, const quadratic_form < long >&A);
  template <> void nugcd < long >(long &x, long &b, long &y, long &a, bool & flag);
  template <> void nucomp_base<long> (quadratic_form<long> & C, const quadratic_form<long> & A, const quadratic_form<long> & B);
  template <> void nudupl_base<long> (quadratic_form<long> & C, const quadratic_form<long> & A);
  template <> long nucomp_base < long >(quadratic_form < long >&C,const quadratic_form < long >&A, const quadratic_form < long >&B, long &G, long &BB);
  template <> long nudupl_base < long > (quadratic_form < long >&C, const quadratic_form < long >&A, long &G, long &B);
  template <> quadratic_form < long >::quadratic_form (const qo_hash_entry < long >&A);
  template <> void quadratic_form < long >::set_current_order (const long &newDelta);
  template <> void quadratic_form < long >::set_current_order (const long &newDelta, const long &newh);
  template <> short quadratic_form<long>::order_signature();
  template <> void quadratic_form < long >::assign_one ();
  template <> void quadratic_form < long >::assign (const qo_hash_entry < long > &B);
  template <> bool quadratic_form < long >::assign_prime (const long &p);
  template <> bool quadratic_form < long >::assign (const long &a2, const long &b2);
  template <> bool quadratic_form < long >::is_equal (const quadratic_form < long >&B) const;
  template <> bool quadratic_form < long >::is_equal (const qo_hash_entry < long >&B) const;
  template <> void nucompcd_base<long> (quadratic_form<long> & C, const quadratic_form<long> & A, const quadratic_form<long> & B);
  template <> void nuduplcd_base<long> (quadratic_form<long> & C, const quadratic_form<long> & A);

  template <> void quadratic_form < long long >::test_form (const char *msg);
  template <> void quadratic_form < long long >::normalize_imag_mult ();
  template <> void quadratic_form < long long >::normalize_imag ();
  template <> void quadratic_form < long long >::normalize_real_mult ();
  template <> void quadratic_form < long long >::normalize_real ();
  template <> void quadratic_form<long long>::normalize_unusual();
  template <> void quadratic_form<long long>::normalize_unusual(quadratic_number<long long> & alpha);
  template <> void quadratic_form < long long >::rho ();
  template <> void quadratic_form < long long >::inverse_rho ();
  template <> void quadratic_form < long long >::rho (quadratic_number<long long> & alpha);
  template <> void quadratic_form < long long >::inverse_rho (quadratic_number<long long> & alpha);
  template <> void quadratic_form < long long >::reduce_imag_mult ();
  template <> void quadratic_form < long long >::reduce_imag_nc ();
  template <> void quadratic_form < long long >::reduce_real ();
  template <> long long multiply_base < long long >(quadratic_form < long long >&C, const quadratic_form < long long >&A, const quadratic_form < long long >&B);
  template <> long long square_base < long long >(quadratic_form < long long >&C, const quadratic_form < long long >&A);
  template <> void nugcd < long long >(long long  &x, long long &b, long long &y, long long &a, bool & flag);
  template <> void nucomp_base<long long> (quadratic_form<long long> & C, const quadratic_form<long long> & A, const quadratic_form<long long> & B);
  template <> void nudupl_base<long long> (quadratic_form<long long> & C, const quadratic_form<long long> & A);
  template <> long long nucomp_base < long long >(quadratic_form < long long >&C,const quadratic_form < long long >&A, const quadratic_form < long long >&B, long long &G, long long &BB);
  template <> long long nudupl_base < long long > (quadratic_form < long long >&C, const quadratic_form < long long >&A, long long &G, long long &B);
  template <> quadratic_form < long long >::quadratic_form (const qo_hash_entry < long long >&A);
  template <> void quadratic_form < long long >::set_current_order (const long long &newDelta);
  template <> void quadratic_form < long long >::set_current_order (const long long &newDelta, const long long &newh);
  template <> short quadratic_form< long long >::order_signature();
  template <> void quadratic_form < long long >::assign_one ();
  template <> void quadratic_form < long long >::assign (const qo_hash_entry < long long > &B);
  template <> bool quadratic_form < long long >::assign_prime (const long long &p);
  template <> bool quadratic_form < long long >::assign (const long long &a2, const long long &b2);
  template <> bool quadratic_form < long long >::is_equal (const quadratic_form < long long >&B) const;
  template <> bool quadratic_form < long long >::is_equal (const qo_hash_entry < long long >&B) const;
  template <> void nucompcd_base<long long> (quadratic_form<long long> & C, const quadratic_form<long long> & A, const quadratic_form<long long> & B);
  template <> void nuduplcd_base<long long> (quadratic_form<long long> & C, const quadratic_form<long long> & A);

} // ANTL

// Unspecialized template definitions.
//#include "../../../src/quadratic/quadratic_form_impl.hpp"

#endif // guard


