/**
 * @file QuadraticOrder.hpp
 * @author Michael Jacobson
 * @brief order in a quadratic number field or hyperelliptic function field
 */

#ifndef ANTL_QUADRATIC_ORDER_H
#define ANTL_QUADRATIC_ORDER_H

#include <vector>

#include <NTL/ZZ.h>
#include <NTL/RR.h>

#include <ANTL/common.hpp>
#include <ANTL/Interface/OrderInvariants.hpp>

// Arithmetic classes
#include <ANTL/Quadratic/Reduce/ReduceStrategy.hpp>
#include <ANTL/Quadratic/Reduce/ReducePlainImag.hpp>
#include <ANTL/Quadratic/Reduce/ReducePlainReal.hpp>
#include <ANTL/Quadratic/Reduce/ReduceFast.hpp>

#include <ANTL/Quadratic/Multiply/MultiplyStrategy.hpp>
#include <ANTL/Quadratic/Multiply/MultiplyPlain.hpp>
#include <ANTL/Quadratic/Multiply/MultiplyNucomp.hpp>

#include <ANTL/Quadratic/Square/SquareStrategy.hpp>
#include <ANTL/Quadratic/Square/SquarePlain.hpp>
#include <ANTL/Quadratic/Square/SquareNudupl.hpp>

#include <ANTL/Quadratic/Cube/CubeStrategy.hpp>
#include <ANTL/Quadratic/Cube/CubePlain.hpp>
#include <ANTL/Quadratic/Cube/CubeNucube.hpp>
#include <ANTL/Quadratic/Cube/CubeMulSqr.hpp>

/*
// Class group classes
#include <ANTL/Quadratic/Invariants/qo_class_group.hpp>
#include <ANTL/Quadratic/Invariants/qo_class_group_bjt.hpp>
#include <ANTL/Quadratic/Invariants/qo_class_group_bs.hpp>
#include <ANTL/Quadratic/Invariants/qo_class_group_subexp.hpp>

// Regulator classes
#include <ANTL/Quadratic/Invariants/qo_regulator.hpp>
#include <ANTL/Quadratic/Invariants/qo_regulator_cfrac.hpp>
#include <ANTL/Quadratic/Invariants/qo_regulator_bsgs.hpp>
#include <ANTL/Quadratic/Invariants/qo_regulator_lenstra.hpp>
#include <ANTL/Quadratic/Invariants/qo_regulator_D16.hpp>
#include <ANTL/Quadratic/Invariants/qo_regulator_subexp.hpp>
*/

namespace ANTL {

using namespace ANTL;

  template < class T > class QuadraticOrder;

  template < class T > class ReduceStrategy;
  template < class T > class ReducePlainImag;
  template < class T > class ReducePlainReal;
  template < class T > class ReduceFast;

  template < class T > class MultiplyStrategy;
  template < class T > class MultiplyPlain;
  template < class T > class MultiplyNucomp;

  template < class T > class SquareStrategy;
  template < class T > class SquarePlain;
  template < class T > class SquareNudupl;

  template < class T > class CubeStrategy;
  template < class T > class CubePlain;
  template < class T > class CubeNucube;
  template < class T > class CubeMulSqr;

  /*
  template < class T >
  QuadraticOrder<T> & randomImaginaryOrder(int size, bool prime);

  template < class T >
  QuadraticOrder<T> & randomUnusualOrder(int size, bool prime);

  template < class T >
  QuadraticOrder<T> & randomRealOrder(int size, bool prime);
  */

  template < class T >
  bool operator == (const QuadraticOrder < T > &QO1, const QuadraticOrder < T > &QO2);

  template < class T >
  bool operator != (const QuadraticOrder < T > &QO1, const QuadraticOrder < T > &QO2);

  template < class T >
  std::istream & operator >> (std::istream & in, QuadraticOrder < T > &QO);

  template < class T >
  std::ostream & operator << (std::ostream & out, const QuadraticOrder < T > &QO);


  /**
   * @brief Order in a quadratic number or function field
   * @remarks This class represents a quadratic order or a hyperelliptic curve,
   * depending on T.  If T is an integer type, then the order
   *    T + (D + sqrt(D))/2 T
   * is represented, where D is the discrinant of the order.  If T represents
   * polynomials over a finite field, then the order
   *    T + rho T
   * is represented, where rho is a root of y^2 + h y + D.
   * Either D or h and D define the order (depending on T).
   *
   * The primary functions of this class are for computing computing invariants
   * of the order related to the ideal class group and infrastructure.  For example,
   * routines for computing the structure and order of the class group, regulator,
   * and fundamental unit are provided, as well as routines for solving
   * the various discrete logarithm problems in this setting.
   *
   * The template parameter DTYPE denotes the type to use for ideal distances in
   * infrastructure computations.  For function fields, this may either be ZZ (default)
   * or long.  For number fields, this can be any of double, quad_float, RR, fpk_rep,
   * compact_representation<T>, or mpfi_t.
   *
   * Currently, the template parameter T has been instantiated for the following types:
   *    base types:
   *       long --- order in quadratic number field (word sized D)
   *    NTL:
   *       ZZ --- order in a quadratic number field (arbitrary sized D)
   *       ZZ_pX --- hyperelliptic function field over Fp (char <> 2)
   *       ZZ_pEX --- hyperelliptic function field over Fq (char <> 2)
   *       zz_pX --- hyperelliptic function field over Fp (char <> 2, p < 2^64)
   *       zz_pEX --- hyperelliptic function field over Fq (char <> 2, p < 2^64)
   *       GF2EX --- hyperelliptic function field (char = 2)
   */

template < class T > class QuadraticOrder : public IOrder<T,NTL::RR> {
  private:

    // TODO: Rename all instances of Delta (and discriminant, etc) to delta throughout codebase
    T hx;    // definition of the order:  y^2 + h(x) y = Delta
    T Delta;
    long g;  // genus (if an order in a function field)

    // Pointers to objects for ideal arithmetic

    ReduceStrategy<T>   *red_best;
    ReducePlainImag<T>  *red_plain_imag;
    ReducePlainReal<T>  *red_plain_real;
    ReduceFast<T>       *red_fast;

    MultiplyStrategy<T> *mul_best;
    MultiplyPlain<T>    *mul_plain;
    MultiplyNucomp<T>   *mul_nucomp;

    SquareStrategy<T>   *sqr_best;
    SquarePlain<T>      *sqr_plain;
    SquareNudupl<T>     *sqr_nudupl;

    CubeStrategy<T>     *cube_best;
    CubePlain<T>        *cube_plain;
    CubeNucube<T>       *cube_nucube;
    CubeMulSqr<T>       *cube_mulsqr;

    //
    // invariants and associated objects objects
    //

    /*
    // Lfunction approximations
    qo_lfunction<T> Lfunc;


    // regulator
    DTYPE R;            // regulator
    bool Rconditional;	// true if correctness of R relies on ERH

    qo_regulator<T,DTYPE> *R_best;
    qo_regulator_cfrac<T,DTYPE> R_cfrac;
    qo_regulator_bsgs<T,DTYPE> R_bsgs;
    qo_regulator_lenstra<T,DTYPE> R_lenstra;
    qo_regulator_D16<T,DTYPE> R_D16;
    qo_regulator_subexp<T,DTYPE> R_subexp;


    // class group
    S h;                          // class number
    vector<S> CL;                 // class group invariants
    vector < qi_class<T> > gens;  // system of generators
    bool hconditional;            // true if correctness of h relies on ERH

    qo_class_group<T,S,DTYPE> *CL_best;
    qo_class_group<T,S,DTYPE> CL_bjt;
    qo_class_group<T,S,DTYPE> CL_bs;
    qo_class_group<T,S,DTYPE> CL_subexp;
    */



  public:
    //
    // constructors and destructor
    //

    QuadraticOrder (const T & D, const T & newh);
    QuadraticOrder (const T & D);
    ~QuadraticOrder ();



    //
    // initialization
    //

    void verbose (long state);
    void verification (long level);


    //
    // random order
    //

    /*
    friend QuadraticOrder<T> & randomImaginaryOrder(long size, bool prime);
    friend QuadraticOrder<T> & randomUnusualOrder(long size, bool prime);
    friend QuadraticOrder<T> & randomRealOrder(long size, bool prime);
    */


    //
    // access functions
    //

    T getH () const
    {
      return hx;
    }

    T getDiscriminant () const
    {
      return Delta;
    }

    long getGenus () const
    {
      return g;
    }

    //
    // comparisons
    //

    bool IsEqual (const QuadraticOrder < T > &QO2) const;

    friend bool operator == < T > (const QuadraticOrder < T > &QO1, const QuadraticOrder < T > &QO2);
    friend bool operator != < T > (const QuadraticOrder < T > &QO1, const QuadraticOrder < T > &QO2);

    //
    // basic functions
    //

    bool IsImaginary () const;
    bool IsUnusual () const;
    bool IsReal () const;

  //Inherited virtual methods


  virtual std::vector<NTL::ZZ> class_group() {
    std::vector<NTL::ZZ> classGroup = {ZZ(1), ZZ(2)};
    return classGroup;
  }

  virtual NTL::ZZ class_number() {
    ZZ classNumber = ZZ(0);
    return classNumber;
  }

  virtual std::vector<T> unit_group() {
    T object1, object2;
    std::vector<T> unitGroup = {object1, object2};
    return unitGroup;
  }

  virtual NTL::RR regulator() {
    NTL::RR regulator;
    return regulator;
  }

    /*
    bool is_R_conditional () const
    {
      return Rconditional;
    };

    bool is_h_conditional () const
    {
      return hconditional;
    };

    bool is_CL_conditional() const
    {
      return hconditional;
    }
    */


    //
    // access to basic ideal arithmetic objects (optimal choice depends on the order)
    //

    void set_red_best(ReduceStrategy<T> &A);
    void set_red_plain_imag(ReducePlainImag<T> &A);
    void set_red_plain_real(ReducePlainReal<T> &A);
    void set_red_fast(ReduceFast<T> &A);

    void set_mul_best(MultiplyStrategy<T> &A);
    void set_mul_plain(MultiplyPlain<T> &A);
    void set_mul_nucomp(MultiplyNucomp<T> &A);

    void set_sqr_best(SquareStrategy<T> &A);
    void set_sqr_plain(SquarePlain<T> &A);
    void set_sqr_nudupl(SquareNudupl<T> &A);

    void set_cube_best(CubeStrategy<T> &A);
    void set_cube_plain(CubePlain<T> &A);
    void set_cube_nucube(CubeNucube<T> &A);
    void set_cube_mulsqr(CubeMulSqr<T> &A);

    ReduceStrategy<T> *   get_red_best();
    ReducePlainImag<T> *  get_red_plain_imag();
    ReducePlainReal<T> *  get_red_plain_real();
    ReduceFast<T> *       get_red_fast();

    MultiplyStrategy<T> * get_mul_best();
    MultiplyPlain<T> *    get_mul_plain();
    MultiplyNucomp<T> *   get_mul_nucomp();

    SquareStrategy<T> *   get_sqr_best();
    SquarePlain<T> *      get_sqr_plain();
    SquareNudupl<T> *     get_sqr_nudupl();

    CubeStrategy<T> *     get_cube_best();
    CubePlain<T> *        get_cube_plain();
    CubeNucube<T> *       get_cube_nucube();
    CubeMulSqr<T> *       get_cube_mulsqr();

    //qo_reduce<T> & reduce_method (long method = 0);
    //qo_multiply<T> & multiply_method (long method = 0);
    //qo_square<T> & square_method (long method = 0);
    //qo_cube<T> & cube_method (long method = 0);



    /*
    //
    // L functions.  Output is templated, should be one of double, RR, or mpfi_t
    //

    void set_unconditional()
    {
      unconditional = true;
    };

    void use_Lfunction_tables(const T & H);
    void set_Lfunction_global(const T & H);


    // precise value via the analytic class number formula, computes regulator
    // and class number first

    template < class RTYPE >
    void Lfunction (RTYPE & outL);

    LTYPE Lfunction();
 template <> bool QuadraticOrder<GF2EX>::IsEqual (const QuadraticOrder<GF2EX> &QO) const;

    //
    // regulator, class number, class group
    //

    void set_regulator (const DTYPE &newR);
    void set_class_number (const S & newh);


    template < class RTYPE >
    void regulator (RTYPE & outR, long method = 0);

    DTYPE regulator(long method = 0);
    //  compact_representation<T> fundamental_unit(long method = 0);

    S class_number ();
    vector<S> class_group (long method = 0);



    // invariants in other floating point types (for improved accuracy, if desired)

    template < class RTYPE >
    void LDfunction (RTYPE & out);

    template < class RTYPE >
    void LLI (RTYPE & out);
qo_multiply_plain
    template < class RTYPE >
    void LLI_D (RTYPE & out);

    template < class RTYPE >
    void ULI (RTYPE & out);

    template < class RTYPE >
    void ULI_D (RTYPE & out);

    long generate_optimal_Q_cfunc ();

    template < class RTYPE >
    void estimate_C (RTYPE & out, long nQ);



    //
    // verification
    //

    bool verify_regulator (bool unconditional = false);
    bool verify_fundamental_unit();
    bool verify_class_group (long level = 0);



    //
    // access to class group statistics
    //
    S & get_CL(long i) {
      S A;
      set(A);
      if (CL.size() > 0 && i >= 0 && i < CL.size())
	A = CL[i];
      return A;
    };

    S exponent ();

    long p_rank (const S & p);

    vector < qi_class<T> > generators ();

    qi_class<T> get_generator(long i) {
      generators();
      qi_class<T> A;
      A.assign_one();
      if (generators.size() > 0 && i >= 0 && i < generators.size())
	A = gens[i];
      return A;
    };

    long get_rank()  { return CL.size(); };

    T get_pmax()  {
      return fact_base[contributors[fact_base.size()-1]].get_a();
    };

    long get_nump()  {
      return fact_base.size();
    };
    */


    //
    // input/output
    //


    void read_from_file (std::istream & in);
    void write_to_file (std::ostream & out) const;


    friend std::istream & operator >> <T> (std::istream & in, QuadraticOrder < T > &QO);
    friend std::ostream & operator << <T> (std::ostream & out, const QuadraticOrder < T > &QO);
  };



  /*
  template <> bool QuadraticOrder<GF2EX,ZZ,ZZ,RR>::assign (const GF2EX & newf);
  template <> bool QuadraticOrder<GF2EX,ZZ,ZZ,RR>::assign (const GF2EX & newf, const GF2EX & newh);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::assign (const QuadraticOrder < GF2EX,ZZ,ZZ,RR > &QO);
  template <> bool QuadraticOrder<GF2EX,ZZ,ZZ,RR>::is_equal (const QuadraticOrder < GF2EX,ZZ,ZZ,RR > &QO) const;
  template <> bool QuadraticOrder<GF2EX,ZZ,ZZ,RR>::is_imaginary () const;
  template <> bool QuadraticOrder<GF2EX,ZZ,ZZ,RR>::is_real () const;
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::randomize_imag (long size, bool prime);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::randomize_unusual (long size, bool prime);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::randomize_real (long size, bool prime);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::use_Lfunction_tables (const GF2EX & H);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::generate_mu_table (long size_start, long size_stop, long num);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::read_from_file (std::istream & in);
  template <> void QuadraticOrder<GF2EX,ZZ,ZZ,RR>::write_to_file (std::ostream & out) const;


  template <> bool QuadraticOrder<long,long,double,double>::assign (const long &D);
  template <> void QuadraticOrder<long,long,double,double>::assign (const QuadraticOrder<long,long,long,double>&QO);
  template <> bool QuadraticOrder<long,long,double,double>::is_imaginary () const;
  template <> bool QuadraticOrder<long,long,double,double>::is_real () const;
  template <> void QuadraticOrder<long,long,double,double>::randomize_imag (long size, bool prime);
  template <> void QuadraticOrder<long,long,double,double>::randomize_unusual (long size, bool prime);
  template <> void QuadraticOrder<long,long,double,double>::randomize_real (long size, bool prime);
  template <> void QuadraticOrder<long,long,double,double>::use_Lfunction_tables (const long & H);
  template <> void QuadraticOrder<long,long,double,double>::set_Lfunction_global (const long & H);
  template <> void QuadraticOrder<long,long,double,double>::Lfunction (double & out);
  template <> void QuadraticOrder<long,long,double,double>::LDfunction (double & out);
  template <> void QuadraticOrder<long,long,double,double>::LLI (double & out);
  template <> void QuadraticOrder<long,long,double,double>::LLI_D (double & out);
  template <> void QuadraticOrder<long,long,double,double>::ULI (double & out);
  template <> void QuadraticOrder<long,long,double,double>::ULI_D (double & out);
  template <> void QuadraticOrder<long,long,double,double>::estimate_C (double & out, long nQ);
  template <> void QuadraticOrder<long,long,double,double>::Lfunction (mpfi_t & out);
  template <> void QuadraticOrder<long,long,double,double>::LDfunction (mpfi_t & out);
  template <> void QuadraticOrder<long,long,double,double>::LLI (mpfi_t & out);
  template <> void QuadraticOrder<long,long,double,double>::LLI_D (mpfi_t & out);
  template <> void QuadraticOrder<long,long,double,double>::ULI (mpfi_t & out);
  template <> void QuadraticOrder<lhpp:141:26ong,long,double,double>::ULI_D (mpfi_t & out);
  template <> void QuadraticOrder<long,long,double,double>::estimate_C (mpfi_t & out, long nQ);
  template <> void QuadraticOrder<long,long,double,double>::write_to_file(std::ostream & out) const;


  template <> bool QuadraticOrder<ZZ,ZZ,RR,RR>::assign (const ZZ & D);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::assign (const QuadraticOrder<ZZ,ZZ,RR,RR> &QO);
  template <> bool QuadraticOrder<ZZ,ZZ,RR,RR>::is_imaginary () const;
  template <> bool QuadraticOrder<ZZ,ZZ,RR,RR>::is_real () const;
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::randomize_imag (long size, bool prime);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::randomize_unusual (long size, bool prime);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::randomize_real (long size, bool prime);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::use_Lfunction_tables (const ZZ & H);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::set_Lfunction_global (const ZZ & H);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::Lfunction (RR & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::LDfunction (RR & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::LLI (RR & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::LLI_D (RR & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::ULI (RR & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::ULI_D (RR & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::estimate_C (RR & out, long nQ);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::Lfunction (mpfi_t & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::LDfunction (mpfi_t & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::LLI (mpfi_t & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::LLI_D (mpfi_t & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::ULI (mpfi_t & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::ULI_D (mpfi_t & out);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::estimate_C (mpfi_t & out, long nQ);
  template <> void QuadraticOrder<ZZ,ZZ,RR,RR>::write_to_file(std::ostream & out) const;


//
// Class quadratic_order is defined as a wrapper for the above base class
// (hides choice of S, DTYPE, and LTYPE from the quadratic_order's base type (T)
//

template<class T>
class QuadraticOrder<T>: public QuadraticOrder<T, ZZ, ZZ, RR> { };

template<>
class QuadraticOrder<ZZ> : public QuadraticOrder<ZZ, ZZ, RR, RR> { };

template<>
class QuadraticOrder<long> : public QuadraticOrder<long, long, double, double> {
};
*/

  template <>      QuadraticOrder<ZZ>::QuadraticOrder (const ZZ & D);
  template <>      QuadraticOrder<long>::QuadraticOrder (const long & D);

  template <> bool QuadraticOrder<GF2EX>::IsEqual (const QuadraticOrder<GF2EX> &QO) const;
  template <> bool QuadraticOrder<GF2EX>::IsImaginary () const;
  template <> bool QuadraticOrder<GF2EX>::IsReal () const;
  //  template <> QuadraticOrder<GF2EX> & randomImaginaryOrder<GF2EX> (long size, bool prime);
  //  template <> QuadraticOrder<GF2EX> & randomUnusualOrder<GF2EX> (long size, bool prime);
  //  template <> QuadraticOrder<GF2EX> & randomRealOrder<GF2EX> (long size, bool prime);

  template <> bool QuadraticOrder<ZZ>::IsImaginary () const;
  template <> bool QuadraticOrder<ZZ>::IsReal () const;
  //  template <> QuadraticOrder<ZZ> & randomImaginaryOrder<ZZ> (long size, bool prime);
  //  template <> QuadraticOrder<ZZ> & randomUnusualOrder<ZZ> (long size, bool prime);
  //  template <> QuadraticOrder<ZZ> & randomRealOrder<ZZ> (long size, bool prime);

  template <> bool QuadraticOrder<long>::IsImaginary () const;
  template <> bool QuadraticOrder<long>::IsReal () const;
  //  template <> QuadraticOrder<long> & randomImaginaryOrder<long> (long size, bool prime);
  //  template <> QuadraticOrder<long> & randomUnusualOrder<long> (long size, bool prime);
  //  template <> QuadraticOrder<long> & randomRealOrder<long> (long size, bool prime);

} // ANTL

// Unspecialized template definitions.
#include "../../../src/Quadratic/QuadraticOrder_impl.hpp"

#endif // guard
