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

template < class T > class QuadraticIdealBase;
template < class T > class QuadraticOrder;
template < class T > class QuadraticNumber;


// declare templated friend functions
template < class T >
void conjugate (QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> &A);

template < class T >
void mul (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

template < class T >
void mul (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

template < class T >
void sqr (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);

template < class T >
void sqr (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);

template < class T >
void cube (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);

template < class T >
void cube (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);

template < class T >
bool operator == (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

template < class T >
bool operator != (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);

template < class T >
std::istream & operator >> (std::istream & in, QuadraticIdealBase<T> &A);

template < class T >
std::ostream & operator << (std::ostream & out, const QuadraticIdealBase<T> &A);



//
// Class: QuadraticIdealBase<T>
//
// This class represents a reduced binary quadratic ideal (a,b,c)
// whose coefficients are of type T.
//

template < class T > class QuadraticIdealBase
{
 protected:

  // ideal coefficients
  T a;
  T b;
  T c;
  ANTL::QuadraticOrder<T> *QO;


public:
  //
  // constructors and destructor
  //

  QuadraticIdealBase (ANTL::QuadraticOrder<T> & inQO);
  ~QuadraticIdealBase ();


  //
  // assignment
  //

  void assign_one ();
  bool assign_prime (const T & p);
  void assign (const QuadraticIdealBase<T> &B);
  void assign(const T & na, const T & nb, const T & nc);
  QuadraticIdealBase<T> &operator = (const QuadraticIdealBase<T> &A);


  // checks whether ideal coeffs are valid (b^2 + bh - ac = Delta)
  void test_ideal (string msg);



  //
  // access functions
  //

  T get_a () const { return a; };
  T get_b () const { return b; };
  T get_c () const { return c; };
  ANTL::QuadraticOrder<T> * get_QO () const { return QO; };


  //
  // arithmetic operations
  //

  friend void conjugate < T > (QuadraticIdealBase<T> & C, const QuadraticIdealBase<T> &A);
  friend void mul < T > (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);
  friend void sqr < T > (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A);
  friend void cube < T > (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A); 

  //
  // arithmetic with relative generator
  //

  friend void mul < T > (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B);
  friend void sqr < T > (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);
  friend void cube < T > (QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A);

  void reduce();
  void reduce(ANTL::QuadraticNumber<T> & gamma);



  //
  // comparisons
  //

  bool IsOne () const;
  bool IsEqual (const QuadraticIdealBase<T> &B) const;

  friend bool operator == < T > (const QuadraticIdealBase<T> &A,
				 const QuadraticIdealBase<T> &B);

  friend bool operator != < T > (const QuadraticIdealBase<T> &A,
				 const QuadraticIdealBase<T> &B);


  //
  // input/output
  //

  friend std::istream & operator >> < T > (std::istream & in,
					   QuadraticIdealBase<T> &A);

  friend std::ostream & operator << < T > (std::ostream & out,
					   const QuadraticIdealBase<T> &A);
};



//
// Declare specialized methods
//

template <> void QuadraticIdealBase<ZZ>::test_ideal(string msg);
template <> void QuadraticIdealBase<ZZ>::assign_one();
template <> bool QuadraticIdealBase<ZZ>::assign_prime (const ZZ & p);

template <> void QuadraticIdealBase<long>::test_ideal(string msg);
template <> void QuadraticIdealBase<long>::assign_one();
template <> bool QuadraticIdealBase<long>::assign_prime (const long & p);

template <> void QuadraticIdealBase<GF2EX>::test_ideal(string msg);
template <> void QuadraticIdealBase<GF2EX>::assign_one();
template <> bool QuadraticIdealBase<GF2EX>::assign_prime (const GF2EX & p);
template <> void conjugate (QuadraticIdealBase<GF2EX> &C, const QuadraticIdealBase<GF2EX> &A);


// Unspecialized template definitions.
#include "../src/Quadratic/QuadraticIdealBase_impl.hpp"

#endif // guard

