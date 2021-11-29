/**
 * @file QuadraticIdealBase.hpp
 * @author Michael Jacobson
 */

#ifndef QUADRATICIDEALBASE_H
#define QUADRATICIDEALBASE_H

#include <string>
#include <ANTL/common.hpp>
#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <ANTL/Quadratic/QuadraticNumber.hpp>

using namespace ANTL;

namespace ANTL {

  template < class T > class QuadraticIdealBase;
  template < class T > class QuadraticOrder;
  template < class T > class QuadraticNumber;

  // declare templated friend functions
  template <class T> void conjugate (QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> &A);

  template <class T> void mul (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

  template <class T> void mul (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

  template <class T> void sqr (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);

  template <class T> void sqr (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);

  template <class T> void cube (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);

  template <class T> void cube (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);

  template <class T> bool operator == (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

  template <class T> bool operator != (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

  template <class T> std::istream & operator >> (std::istream & in, QuadraticIdealBase<T> &A);

  template <class T> std::ostream & operator << (std::ostream & out, const QuadraticIdealBase<T> &A);

  // Class: QuadraticIdealBase<T>
  //
  // This class represents a primtive integral ideal as a binary quadratic form with coefficients a,b,c of type T.
  template <class T> class QuadraticIdealBase {
    protected:
      // BQF coefficients
      T a;
      T b;
      T c;
      QuadraticOrder<T> *QO;

    public:
      // Constructor(s) and destructor
      QuadraticIdealBase (ANTL::QuadraticOrder<T> & inQO);
      ~QuadraticIdealBase ();

      // Assignment
      void assign_one   ();
      bool assign_prime (const T & p);
      void assign       (const QuadraticIdealBase<T> &B);
      void assign       (const T & na, const T & nb, const T & nc);

      QuadraticIdealBase<T> &operator = (const QuadraticIdealBase<T> &A);

      // Checks whether ideal coeffs are valid (b^2 + bh - ac = Delta)
      void ensure_valid (std::string msg);

      // Getters and Setters
      T get_a () const;
      T get_b () const;
      T get_c () const;
      QuadraticOrder<T> * get_QO () const;

      void set_a  (T x);
      void set_b  (T x);
      void set_c  (T x);
      void set_QO (QuadraticOrder<T> *qo);

      // arithmetic operations
      friend void conjugate < T > (QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> &A);
      friend void mul < T > (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);
      friend void sqr < T > (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);
      friend void cube < T > (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);

      // arithmetic with relative generator
      //   friend void mul < T > (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);
      //   friend void sqr < T > (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);
      //   friend void cube < T > (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);

      void reduce();
      void reduce(ANTL::QuadraticNumber<T> & gamma);

      // comparisons
      bool IsOne () const;
      bool IsEqual (const QuadraticIdealBase<T> &B) const;

      bool is_normal();
      bool is_reduced();

      void normalize();

      friend bool operator == < T > (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);
      friend bool operator != < T > (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

      // input/output
      friend std::istream & operator >> < T > (std::istream & in, QuadraticIdealBase<T> &A);
      friend std::ostream & operator << < T > (std::ostream & out, const QuadraticIdealBase<T> &A);
    };

  // Declare specialized methods
  template <> void QuadraticIdealBase<ZZ>::ensure_valid(std::string msg);
  template <> void QuadraticIdealBase<ZZ>::assign_one();
  template <> bool QuadraticIdealBase<ZZ>::assign_prime (const ZZ & p);
  template <> bool QuadraticIdealBase<ZZ>::is_normal();
  template <> bool QuadraticIdealBase<ZZ>::is_reduced();


  template <> void QuadraticIdealBase<long>::ensure_valid(std::string msg);
  template <> void QuadraticIdealBase<long>::assign_one();
  template <> bool QuadraticIdealBase<long>::assign_prime (const long & p);
  template <> bool QuadraticIdealBase<long>::is_normal();
  template <> bool QuadraticIdealBase<long>::is_reduced();
  template <> void QuadraticIdealBase<ZZ>::normalize();

  template <> void QuadraticIdealBase<GF2EX>::ensure_valid(std::string msg);
  template <> void QuadraticIdealBase<GF2EX>::assign_one();
  template <> bool QuadraticIdealBase<GF2EX>::assign_prime (const GF2EX & p);
  template <> void conjugate (QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A);

} // ANTL

// Unspecialized template definitions.
#include "../src/Quadratic/QuadraticIdealBase_impl.hpp"

#endif // guard

