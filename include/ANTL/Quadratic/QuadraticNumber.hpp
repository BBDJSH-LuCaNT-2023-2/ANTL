/**
 * @file QuadraticNumber.hpp
 * @author Michael Jacobson
 * @brief element of an order in a quadratic number field or hyperelliptic function field
 */

#ifndef ANTL_QUADRATIC_NUMBER_H
#define ANTL_QUADRATIC_NUMBER_H

#include <ANTL/common.hpp>
#include <ANTL/Arithmetic/QQ.hpp>
#include <ANTL/Quadratic/QuadraticOrder.hpp>

using namespace ANTL;

namespace ANTL {
  // Forward Declarations
  template <class T> class QuadraticNumber;

  template <class T> class QuadraticOrder;

  template<class T> void
    clear (QuadraticNumber<T> & x);

  template<class T>
    void
    set (QuadraticNumber<T> & x);

  template<class T>
    void
    assign (QuadraticNumber<T> & z, const QuadraticNumber<T> & x);

  template<class T>
    bool
    IsZero (const QuadraticNumber<T>& x);

  template<class T>
    bool
    IsOne (const QuadraticNumber<T>& x);

  template<class T>
    bool
    IsEqual (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    bool
    IsEqual (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    bool
    IsEqual (const T & n, const QuadraticNumber<T> & x);

  template<class T>
    bool
    IsEqual (const QuadraticNumber<T> & x, const QQ<T> & n);

  template<class T>
    bool
    IsEqual (const QQ<T> & n, const QuadraticNumber<T> & x);

  template<class T>
    bool
    operator == (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    bool
    operator == (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    bool
    operator == (const T & n, const QuadraticNumber<T> & x);

  template<class T>
    bool
    operator == (const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    bool
    operator == (const QQ<T> & q, const QuadraticNumber<T> & x);

  template<class T>
    bool
    operator != (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    bool
    operator != (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    bool
    operator != (const T & n, const QuadraticNumber<T> & x);

  template<class T>
    bool
    operator != (const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    bool
    operator != (const QQ<T> & q, const QuadraticNumber<T> & x);

  template<class T>
    void
    conjugate (QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    void
    inv (QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    void
    add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	 const QuadraticNumber<T> & y);

  template<class T>
    void
    add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n);

  template<class T>
    void
    add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    void
    sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	 const QuadraticNumber<T> & y);

  template<class T>
    void
    sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n);

  template<class T>
    void
    sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    void
    mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	 const QuadraticNumber<T> & y);

  template<class T>
    void
    mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n);

  template<class T>
    void
    mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    void
    div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	 const QuadraticNumber<T> & y);

  template<class T>
    void
    div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n);

  template<class T>
    void
    div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    void
    sqr (QuadraticNumber<T> & z, const QuadraticNumber<T> & x);

  template<class T>
    void
    cube (QuadraticNumber<T> & z, const QuadraticNumber<T> & x);

  template<class T>
    QuadraticNumber<T>
    operator - (const QuadraticNumber<T> & x);

  template<class T>
    QuadraticNumber<T>
    operator + (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    QuadraticNumber<T>
    operator + (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    QuadraticNumber<T>
    operator + (const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    QuadraticNumber<T>
    operator - (const QuadraticNumber<T> & x, const QuadraticNumber<T> & b);

  template<class T>
    QuadraticNumber<T>
    operator - (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    QuadraticNumber<T>
    operator - (const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    QuadraticNumber<T>
    operator * (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    QuadraticNumber<T>
    operator * (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    QuadraticNumber<T>
    operator * (const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    QuadraticNumber<T>
    operator / (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y);

  template<class T>
    QuadraticNumber<T>
    operator / (const QuadraticNumber<T> & x, const T & n);

  template<class T>
    QuadraticNumber<T>
    operator / (const QuadraticNumber<T> & x, const QQ<T> & q);

  template<class T>
    std::istream &
    operator>> (std::istream & in, QuadraticNumber<T> & x);

  template<class T>
    std::ostream &
    operator<< (std::ostream & out, const QuadraticNumber<T> & x);



  /**
   * @brief element of an order in a quadratic number field or hyperelliptic function field
   * @remarks A quadratic number is an element of an order in Q(rho) where rho is a root of an irreducible quadratic
   * polynomial with coefficients of type T is represented as
   *        (a + b \rho)/d
   * where a,b,d \in T.
   *
   * Instances of QuadraticNumber are all normalized to satisfy the following properties:
   *  - gcd(a,b,d) == 1.
   *  - if T is a numerical type (i.e. the element belongs to an order in a quadratic number field), then
   *    d is positive.
   *  - If T is a polynomial type (i.e. the element belongs to an order in a hyperelliptic function field), then
   *    - d is monic
   *    - a is monic when deg(a) > deg(b)+g, and b is monic when deg(a) <= deg(b)+g
   *    where g is the genus of the function field.
   *
   * A QuadraticNumber must be initialized with an existing QuadraticOrder, indicating which order it belongs to.
   * This order is fixed throughout the lifetime of the instance.
   *
   * Currently, the following template parameters T are instantiated:
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
  template<class T>
    class QuadraticNumber
    {
    private:
      QuadraticOrder<T> *QO; /**< order to which the QuadraticNumber belongs */
      T a, b, d; /**< coefficients of the QuadraticNumber */

      static T newA,newB,newD,temp;	// temporary variables for arithmetic operations

      void
      normalize ()
      {
	// remove common factors
	T g = GCD(GCD(a,b),d);
	if (!::IsOne(g)) {
	    ::div(a,a,g);
	    ::div(b,b,g);
	    ::div(d,d,g);
	}

	// normalize leading coefficients
	MakeMonic(d);
	if (deg(a) > (deg(b) + QO->getGenus()))
	  MakeMonic(a);
	else
	  MakeMonic(b);
      }


    public:
      /**
       * @brief Initializes a zero QuadraticNumber in the given QuadraticOrder
       * @param[in] inQO the QuadraticOrder to which this instance belongs
       */
      QuadraticNumber (QuadraticOrder<T> & inQO)
      {
	QO = &inQO;
	::clear (a);
	::clear (b);
	::set (d);
      }

      /**
       * @brief Copy constructor
       * @param[in] n the constant to copy
       */
      QuadraticNumber (QuadraticOrder<T> & inQO, T & n)
      {
	QO = &inQO;
	a = n;
	::clear(b);
	::set(d);
      }

      /**
       * @brief Copy constructor
       * @param[in] q the rational number or function to copy
       */
      QuadraticNumber (QuadraticOrder<T> & inQO, QQ<T> & q)
      {
	QO = &inQO;
	a = q.getNumerator();
	::clear(b);
	d = q.getDenominator();
      }

       /**
       * @brief Copy constructor
       * @param[in] x the QuadraticNumber to copy
       */
      QuadraticNumber (QuadraticNumber<T> & x)
      {
	QO = x.QO;
	a = x.a;
	b = x.b;
	d = x.d;
      }

      ~QuadraticNumber ()
      {
      }



      /**
       * Accessor methods
       */

      const T &
      get_a () const
      {
	return a;
      }

      const T &
      get_b () const
      {
	return b;
      }

      const T &
      get_d () const
      {
	return d;
      }

      QuadraticOrder<T> *
      getQO () const
      {
	return QO;
      }

      /**
       * @brief Computes the norm of this QuadraticNumber
       * @return the norm as a QQ (rational number or function)
       */
      QQ<T> &
      getNorm () const
      {
	QQ<T> q;

	::sqr(newA,a);
	::sqr(temp,b);
	::mul(temp,temp,QO->getDiscriminant());
	::sub(newA,newA,temp);

	::sqr(newD,d);

	q.assign(newA,newD);
	return q;
      }

      /**
       * @brief Computes the trace of this QuadraticNumber
       * @return the trace as a QQ (rational number or function)
       */
      QQ<T> &
      getTrace () const
      {
	QQ<T> q;

	::add(newA,a,a);
	q.assign(newA,d);
	return q;
      }


      /**
       *  Assignments and mutator methods.  Note that the QO value may not be changed.
       **/

      void
      set_a (const T & inA)
      {
	a = inA;
	normalize ();
      }

      void
      set_b (const T & inB)
      {
	b = inB;
	normalize ();
      }

      void
      set_d (const T & inD)
      {
	d = inD;
	normalize ();
      }

      /**
       * @brief Sets z to zero
       * @param[out] z QuadraticNumber to be set to zero
       */
      friend void
      clear (QuadraticNumber<T> & z)
      {
	::clear (z.a);
	::clear (z.b);
	::set (z.d);
      }

      /**
       * @brief Sets z to one
       * @param[out] z QuadraticNumber to be set to one
       */
      friend void
      set (QuadraticNumber<T> & z)
      {
	::set (z.a);
	::clear (z.b);
	::set (z.d);
      }

      /**
       * @brief Sets this QuadraticNumber equal to x
       * @param[in] x value to give the QuadraticNumber
       * @pre this QuadraticNumber and x must belong to the same QuadraticOrder
       */
      void
      assign (const QuadraticNumber<T> & x)
      {
	if (QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }
	a = x.a;
	b = x.b;
	d = x.d;
      }

      /**
       * @brief Sets this QuadraticNumber equal to n
       * @param[in] n value to give the QuadraticNumber
       */
      void
      assign (const T & n)
      {
	a = n;
	::clear (b);
	::set (d);
      }

      /**
       * @brief Sets this QuadraticNumber equal to q
       * @param[in] q value to give the QuadraticNumber
       */
      void
      assign (const QQ<T> & q)
      {
	a = q.getNumerator ();
	::clear (b);
	d = q.getDenominator ();
      }

      /**
       * @brief Sets the QuadraticNumber z equal to x
       * @param[out] z result of assignment
       * @param[in] x value to be assigned
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      assign (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
	if (z.QO != x.QO)
	{
	    // TODO: THROW AN EXCEPTION HERE!!!
	}

	z.a = x.a;
	z.b = x.b;
	z.d = x.d;
      }

      /**
       * @brief Assignment operator
       * @param[in] x value to give the QuadraticNumber
       */
      const QuadraticNumber<T> &
      operator= (const QuadraticNumber<T> & x)
      {
	this->assign (x);
	return *this;

      }



      /**
       * Comparison methods
       **/

      /**
       * @brief Tests whether this QuadraticNumber is equal to zero
       * @return True if the QuadraticNumber is zero
       */
      bool
      isZero () const
      {
	return (::IsZero (a) && ::IsZero (b) && ::IsOne (d));
      }
      ;

      /**
       * @brief Tests whether this QuadraticNumber is equal to one
       * @return True if the QuadraticNumber is one
       */
      bool
      isOne () const
      {
	return (::IsOne (a) && ::IsZero (b) && ::IsOne (d));
      }

      /**
       * @brief Tests whether x is equal to zero
       * @param[in] x QuadraticNumber to test
       * @return True if x is zero
       */
      friend bool
      IsZero (const QuadraticNumber<T> & x)
      {
	return (::IsZero (x.a) && ::IsZero (x.b) && ::IsOne (x.d));
      }

      /**
       * @brief Tests whether x is equal to one
       * @param[in] x QuadraticNumber to test
       * @return True if x is one
       */
      friend bool
      IsOne (const QuadraticNumber<T> & x)
      {
	return (::IsOne (x.a) && ::IsZero (x.b) && ::IsOne (x.d));
      }

      /**
       * @brief Determines whether this QuadraticNumber is integral (norm is an integer)
       * @return true if this is integral
       */
      bool
      isIntegral () const
      {
	return IsInteger(getNorm());
      }

      /**
       * @brief Determines whether this QuadraticNumber is a unit
       * @return true if this is a unit (norm = 1)
       */
      bool
      isUnit () const
      {
	// Note that in the function field case, the normalization we have imposed implies that if this is
	// a unit, then the norm will be equal to 1 (not just a constant)

	return deg(getNorm()) == 0;
      }

      /**
       * @brief Tests whether this QuadraticNumber is equal to y
       * @param[in] y QuadraticNumber to compare
       * @return True if the QuadraticNumber is equal to y
       */
      bool
      isEqual (const QuadraticNumber<T> & x)
      {
	return (QO == x.QO && a == x.a && b = x.b && d = x.d);
      }

      /**
       * @brief Tests whether this QuadraticNumber is equal to n
       * @param[in] n constant to compare
       * @return True if the QuadraticNumber is equal to n
       */
      bool
      isEqual (const T & n)
      {
	return (a == n && IsZero (b) && IsOne (d));
      }

      /**
       * @brief Tests whether this QuadraticNumber is equal to q
       * @param[in] q rational number or function to compare
       * @return True if the QuadraticNumber is equal to q
       */
      bool
      isEqual (const QQ<T> & q)
      {
	return (a == q.getNumerator () && IsZero (b) && d == q.getDenominator ());
      }

      /**
       * @brief Tests whether x is equal to y
       * @param[in] x QuadraticNumber to compare
       * @param[in] y QuadraticNumber to compare
       * @return True if x is equal to y
       */
      friend bool
      IsEqual (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
	return (x.QO == y.QO && x.a == y.a && x.b = y.b && x.d = y.d);
      }

      /**
       * @brief Tests whether x is equal to x
       * @param[in] x QuadraticNumber to compare
       * @param[in] n constant to compare
       * @return True if x is equal to n
       */
      friend bool
      IsEqual (const QuadraticNumber<T> & x, const T & n)
      {
	return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

      /**
       * @brief Tests whether x is equal to n
       * @param[in] n constant to compare
       * @param[in] x QuadraticNumber to compare
       * @return True if x is equal to n
       */
      friend bool
      IsEqual (const T & n, const QuadraticNumber<T> & x)
      {
	return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

      /**
       * @brief Tests whether x is equal to q
       * @param[in] x QuadraticNumber to compare
       * @param[in] q rational number or function to compare
       * @return True if x is equal to q
       */
      friend bool
      IsEqual (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

      /**
       * @brief Tests whether x is equal to q
       * @param[in] q rational number or function to compare
       * @param[in] x QuadraticNumber to compare
       * @return True if x is equal to q
       */
      friend bool
      IsEqual (const QQ<T> & q, const QuadraticNumber<T> & x)
      {
	return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

      /**
       * @brief Equality test operator
       * @param[in] x QuadraticNumber to compare
       * @param[in] y QuadraticNumber to compare
       * @return True if x is equal to y
       */
      friend bool
      operator == (const QuadraticNumber<T> & x,
		      const QuadraticNumber<T> & y)
      {
	return (x.QO == y.QO && x.a == y.a && x.b = y.b && x.d = y.d);
      }

      /**
       * @brief Equality test operator
       * @param[in] x QuadraticNumber to compare
       * @param[in] n constant to compare
       * @return True if x is equal to n
       */
      friend bool
      operator == (const QuadraticNumber<T> & x, const T & n)
      {
	return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

      /**
       * @brief Equality test operator
       * @param[in] n constant to compare
       * @param[in] x QuadraticNumber to compare
       * @return True if x is equal to n
       */
      friend bool
      operator == (const T & n, const QuadraticNumber<T> & x)
      {
	return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

      /**
       * @brief Equality test operator
       * @param[in] x QuadraticNumber to compare
       * @param[in] q rational number or function to compare
       * @return True if x is equal to q
       */
      friend bool
      operator == (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

      /**
       * @brief Equality test operator
       * @param[in] q rational number or function to compare
       * @param[in] x QuadraticNumber to compare
       * @return True if x is equal to q
       */
      friend bool
      operator == (const QQ<T> & q, const QuadraticNumber<T> & x)
      {
	return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

      /**
       * @brief Inequality test operator
       * @param[in] x QuadraticNumber to compare
       * @param[in] y QuadraticNumber to compare
       * @return True if x is not equal to y
       */
      friend bool
      operator != (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
	return !(x == y);
      }

      /**
       * @brief Inequality test operator
       * @param[in] x QuadraticNumber to compare
       * @param[in] n constant to compare
       * @return True if x is not equal to n
       */
      friend bool
      operator != (const QuadraticNumber<T> & x, const T & n)
      {
	return !(x == n);
      }

      /**
       * @brief Inequality test operator
       * @param[in] n constant to compare
       * @param[in] x QuadraticNumber to compare
       * @return True if x is not equal to n
       */
      friend bool
      operator != (const T & n, const QuadraticNumber<T> & x)
      {
	return !(x == n);
      }

      /**
       * @brief Inequality test operator
       * @param[in] x QuadraticNumber to compare
       * @param[in] q rational number or function to compare
       * @return True if x is not equal to q
       */
      friend bool
      operator != (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	return !(x == q);
      }

      /**
       * @brief Inequality test operator
       * @param[in] q rational number or function to compare
       * @param[in] x QuadraticNumber to compare
       * @return True if x is not equal to q
       */
      friend bool
      operator != (const QQ<T> & q, const QuadraticNumber<T> & x)
      {
	return !(x == q);
      }



      /**
       * Arithmetic operations
       **/

      /**
       * @brief Negates the QuadraticNumber
       */
      void
      negate ()
      {
	a = -a;
	b = -b;
	normalize();
      }

      /**
       * @brief Inverts the QuadraticNumber
       */
      void
      invert ()
      {
	// ((a + b rho) / d)^-1 = (ad - bd rho) / (a^2 - b^2 Delta)
	::mul(newA,a,d);

	::mul(newB,b,d);

	::sqr(newD,a);
	::sqr(temp,b);
	::mul(temp,temp,QO->getDiscriminant());
	::sub(newD,newD,temp);

	a = newA;
	b = -newB;
	d = newD;

	normalize ();
      }



      /**
       * @brief Sets the QuadraticNumber z equal to the conjugate of x
       * @param[out] z conjugate of x
       * @param[in] x value to be conjugated
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      conjugate (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
	if (z.QO != x.QO) {
	  // TODO:  THROW AN EXCEPTION!!!
	}

	z.a = x.a;
	z.b = -x.b;
	z.d = x.d;
      }

      /**
       * @brief Sets the QuadraticNumber z equal to the inverse of x
       * @param[out] z inverse of x
       * @param[in] x value to be inverted
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      inv (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
	if (z.QO != x.QO) {
	  // TODO:  THROW AN EXCEPTION!!!
	}

	// ((a + b rho) / d)^-1 = (ad - bd rho) / (a^2 - b^2 Delta)
	::mul(z.a,x.a,x.d);

	::mul(z.b,x.b,x.d);

	::sqr(z.d,x.a);
	::sqr(temp,x.b);
	::mul(temp,temp,x.QO->getDiscriminant());
	::sub(z.d,z.d,temp);

	z.normalize ();
      }

      /**
       * @brief Computes the sum of x and y
       * @param[out] z = x + y
       * @param[in] x first summand
       * @param[in] y second summand
       * @pre z, x, and y must belong to the same QuadraticOrder
       */
      friend void
      add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QuadraticNumber<T> & y)
      {
	if (z.QO != x.QO || z.QO != y.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	if (x.isZero())
	  z.assign (y);
	else if (y.isZero())
	  z.assign (x);
	else
	  {
	    //     y.d * x.a + x.d * y.a + rho * (y.d * x.b + x.d * y.b)
	    // z = ---------------------------------------------------------
	    //                     x.d * y.d

	    ::mul(newA,x.a,y.d);
	    ::mul(temp,y.a,x.d);
	    ::add(newA,newA,temp);

	    ::mul(newB,x.b,y.d);
	    ::mul(temp,y.b,x.d);
	    ::add(newB,newB,temp);

	    ::mul(newD,x.d,y.d);

	    z.a = newA;
	    z.b = newB;
	    z.d = newD;
	    z.normalize();
	  }
      }

      /**
       * @brief Computes the sum of x and n (a constant)
       * @param[out] z = x + n
       * @param[in] x first summand
       * @param[in] n second summand (a constant)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	if (x.isZero())
	  z.assign (n);
	else if (IsZero(n))
	  z.assign (x);
	else
	  {
	    //     (x.a + x.d * n) + x.b * rho
	    // z = ---------------------------------------------------------
	    //                x.d

	    ::mul(temp,x.d,n);
	    ::add(z.a,z.a,temp);
	    z.b = x.b;
	    z.d = x.d;
	    z.normalize();
	  }

      }

      /**
       * @brief Computes the sum of x and q (a rational number or function)
       * @param[out] z = x + q
       * @param[in] x first summand
       * @param[in] q second summand (a rational number or function)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QQ<T> & q)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	if (x.isZero())
	  z.assign (q);
	else if (IsZero(q))
	  z.assign (x);
	else
	  {
	    //     q.d * x.a + x.d * q.n + rho * (q.d * x.b)
	    // z = ---------------------------------------------------------
	    //                     x.d * q.d

	    ::mul(newA,x.a,q.getDenominator());
	    ::mul(temp,q.getNumerator(),x.d);
	    ::add(newA,newA,temp);

	    ::mul(newB,x.b,q.getDenominator());

	    ::mul(newD,x.d,q.getDenominator());

	    z.a = newA;
	    z.b = newB;
	    z.d = newD;
	    z.normalize();
	  }
      }


      /**
       * @brief Computes the difference of x and y
       * @param[out] z = x - y
       * @param[in] x first term
       * @param[in] y second term
       * @pre z, x, and y must belong to the same QuadraticOrder
       */
      friend void
      sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QuadraticNumber<T> & y)
      {
	if (z.QO != x.QO || z.QO != y.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	if (x.isZero()) {
	    z.assign(y);
	    z.negate();
	}
	else if (y.isZero())
	  z.assign (x);
	else
	  {
	    //     y.d * x.a - x.d * y.a + rho * (y.d * x.b - x.d * y.b)
	    // z = ---------------------------------------------------------
	    //                     x.d * y.d

	    ::mul(newA,x.a,y.d);
	    ::mul(temp,y.a,x.d);
	    ::sub(newA,newA,temp);

	    ::mul(newB,x.b,y.d);
	    ::mul(temp,y.b,x.d);
	    ::sub(newB,newB,temp);

	    ::mul(newD,x.d,y.d);

	    z.a = newA;
	    z.b = newB;
	    z.d = newD;
	    z.normalize();
	  }
      }

      /**
       * @brief Computes the difference of x and n (a constant)
       * @param[out] z = x - n
       * @param[in] x first term
       * @param[in] n second term (a constant)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	if (x.isZero()) {
	  z.assign(n);
	  z.negate();
	}
	else if (IsZero(n))
	  z.assign (x);
	else
	  {
	    //     (x.a - x.d * n) + x.b * rho
	    // z = ---------------------------------------------------------
	    //                x.d

	    ::mul(temp,x.d,n);
	    ::sub(z.a,z.a,temp);
	    z.b = x.b;
	    z.d = x.d;
	    z.normalize();
	  }

      }

      /**
       * @brief Computes the difference of x and q (a rational number or function)
       * @param[out] z = x - q
       * @param[in] x first term
       * @param[in] q second term (a rational number or function)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QQ<T> & q)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	if (x.isZero()) {
	  z.assign (q);
	  z.negate();
	}
	else if (IsZero(q))
	  z.assign (x);
	else
	  {
	    //     q.d * x.a - x.d * q.n + rho * (q.d * x.b)
	    // z = ---------------------------------------------------------
	    //                     x.d * q.d

	    ::mul(newA,x.a,q.getDenominator());
	    ::mul(temp,q.getNumerator(),x.d);
	    ::sub(newA,newA,temp);

	    ::mul(newB,x.b,q.getDenominator());

	    ::mul(newD,x.d,q.getDenominator());

	    z.a = newA;
	    z.b = newB;
	    z.d = newD;
	    z.normalize();
	  }
      }


      /**
       * @brief Computes the product of x and y
       * @param[out] z = x * y
       * @param[in] x first term
       * @param[in] y second term
       * @pre z, x, and y must belong to the same QuadraticOrder
       */
      friend void
      mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QuadraticNumber<T> & y)
      {
	if (z.QO != x.QO || z.QO != y.QO) {
	  // TODO:  THROW AN EXCEPTION
        }

	// z = (x.a y.a + x.b y.b N(rho)) + rho(x.a y.b + y.a x.b + x.b y.b Tr(rho))
	//     ---------------------------------------------------------------------
	//                              x.d y.d

	::mul(newA,x.a,y.a);
	::mul(temp,x.b,y.b);
	::mul(temp,temp,z.QO->getDiscriminant());
	::add(newA,newA,temp);

	::mul(newB,x.a,y.b);
	::mul(temp,x.b,y.a);
	::add(newB,newB,temp);

	::mul(newD,x.d,y.d);

	z.a = newA;
	z.b = newB;
	z.d = newD;
	z.normalize ();
      }

      /**
       * @brief Computes the product of x and n (a constant)
       * @param[out] z = x * n
       * @param[in] x first term
       * @param[in] n second term (a constant)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	::mul(z.a,x.a,n);
	::mul(z.b,x.b,n);
	z.d = x.d;
	z.normalize();
      }

      /**
       * @brief Computes the product of x and q (a rational number or function)
       * @param[out] z = x * q
       * @param[in] x first term
       * @param[in] q second term (a rational number or function)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QQ<T> & q)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	::mul(z.a,x.a,q.getNumerator());
	::mul(z.b,x.b,q.getNumerator());
	::mul(z.d,z.d,q.getDenominator());
	z.normalize();
      }

      /**
       * @brief Computes the quotient of x and y
       * @param[out] z = x / y
       * @param[in] x first term
       * @param[in] y second term
       * @pre z, x, and y must belong to the same QuadraticOrder
       */
      friend void
      div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QuadraticNumber<T> & y)
      {
	if (z.QO != x.QO || z.QO != y.QO) {
	  // TODO:  THROW AN EXCEPTION
        }

	// z = y.d (x.a y.a - x.b y.b N(rho)) - y.d (x.b y.a - x.a y.b) rho
	//     ------------------------------------------------------------
	//                   x.d (y.a^2 - y.b^2 N(rho)

	::mul(newA,x.a,y.a);
	::mul(temp,x.b,y.b);
	::mul(temp,temp,z.QO->getDiscriminant());
	::sub(newA,newA,temp);
	::mul(newA,newA,y.d);

	::mul(newB,x.b,y.a);
	::mul(temp,x.a,y.b);
	::sub(newB,newB,temp);
	::mul(newB,newB,y.d);

	::sqr(newD,y.a);
	::sqr(temp,y.b);
	::mul(temp,temp,z.QO->getDiscriminant());
	::mul(newD,temp,x.d);

	z.a = newA;
	z.b = newB;
	z.d = newD;
	z.normalize ();
      }

      /**
       * @brief Computes the quotient of x and n (a constant)
       * @param[out] z = x / n
       * @param[in] x first term
       * @param[in] n second term (a constant)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	z.a = x.a;
	z.b = x.b;
	::mul(z.d,x.d,n);
	z.normalize();
      }

      /**
       * @brief Computes the quotient of x and q (a rational number or function)
       * @param[out] z = x / q
       * @param[in] x first term
       * @param[in] q second term (a rational number or function)
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x,
	   const QQ<T> & q)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	::mul(z.a,x.a,q.getDenominator());
	::mul(z.b,x.b,q.getDenominator());
	::mul(z.d,x.d,q.getNumerator());
	z.QO = x.QO;
	z.normalize();
      }

      /**
       * @brief Computes the square of x
       * @param[out] z = x^2
       * @param[in] x value to square
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      sqr (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	// z = (x.a^2 + y.b^2 N(rho)) + rho( 2 x.a y.b + x.b y.b Tr(rho))
	//     ----------------------------------------------------------
	//                         x.d^2

	::sqr(newA,x.a);
	::sqr(temp,x.b);
	::mul(temp,temp,z.QO->getDiscriminant());
	::add(newA,newA,temp);

	::mul(newB,x.a,x.b);
	::add(newB,newB,newB);

	::sqr(newD,x.d);

	z.a = newA;
	z.b = newB;
	z.d = newD;
	z.normalize ();
      }

      /**
       * @brief Computes the cube of x
       * @param[out] z = x^3
       * @param[in] x value to cube
       * @pre z and x must belong to the same QuadraticOrder
       */
      friend void
      cube (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
	if (z.QO != x.QO)
	  {
	    // TODO:  THROW AN EXCEPTION!!!
	  }

	QuadraticNumber<T> zz;
	sqr(zz,x);
	mul(z,zz,x);
      }


      friend QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & a)
      {
	QuadraticNumber<T> c (a);
	c.negate ();
	return c;
      }

      friend QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
	QuadraticNumber<T> c;
	add (c, x, y);
	return c;
      }

      friend QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const T & n)
      {
	QuadraticNumber<T> c;
	add (c, x, n);
	return c;
      }

      friend QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	QuadraticNumber<T> c;
	add (c, x, q);
	return c;
      }

      friend QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
	QuadraticNumber<T> c;
	sub (c, x, y);
	return c;
      }

      friend QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const T & n)
      {
	QuadraticNumber<T> c;
	sub (c, x, n);
	return c;
      }

      friend QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	QuadraticNumber<T> c;
	sub (c, x, q);
	return c;
      }

      friend QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
	QuadraticNumber<T> c;
	mul (c, x, y);
	return c;
      }

      friend QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const T & n)
      {
	QuadraticNumber<T> c;
	mul (c, x, n);
	return c;
      }

      friend QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	QuadraticNumber<T> c;
	mul (c, x, q);
	return c;
      }

      friend QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
	QuadraticNumber<T> c;
	div (c, x, y);
	return c;
      }

      friend QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & a, const T & n)
      {
	QuadraticNumber<T> c;
	div (c, a, n);
	return c;
      }

      friend QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
	QuadraticNumber<T> c;
	div (c, x, q);
	return c;
      }

      const QuadraticNumber<T> &
      operator += (const QuadraticNumber<T> & x)
      {
	add (*this, *this, x);
	return *this;
      }

      const QuadraticNumber<T> &
      operator += (const T & n)
      {
	add (*this, *this, n);
	return *this;
      }

      const QuadraticNumber<T> &
      operator += (const QQ<T> & q)
      {
	add (*this, *this, q);
	return *this;
      }

      const QuadraticNumber<T> &
      operator -= (const QuadraticNumber<T> & x)
      {
	sub (*this, *this, x);
	return *this;
      }

      const QuadraticNumber<T> &
      operator -= (const T & n)
      {
	sub (*this, *this, n);
	return *this;
      }

      const QuadraticNumber<T> &
      operator -= (const QQ<T> & q)
      {
	sub (*this, *this, q);
	return *this;
      }

      const QuadraticNumber<T> &
      operator *= (const QuadraticNumber<T> & x)
      {
	mul (*this, *this, x);
	return *this;
      }

      const QuadraticNumber<T> &
      operator *= (const T & n)
      {
	mul (*this, *this, n);
	return *this;
      }

      const QuadraticNumber<T> &
      operator *= (const QQ<T> & q)
      {
	mul (*this, *this, q);
	return *this;
      }

      const QuadraticNumber<T> &
      operator /= (const QuadraticNumber<T> & x)
      {
	div (*this, *this, x);
	return *this;
      }

      const QuadraticNumber<T> &
      operator /= (const T & n)
      {
	div (*this, *this, n);
	return *this;
      }

      const QuadraticNumber<T> &
      operator /= (const QQ<T> & q)
      {
	div (*this, *this, q);
	return *this;
      }



      /**
       *  input / output
       **/

      friend std::istream &
      operator>> (std::istream & in, QuadraticNumber<T> & x)
      {
	    char c;

	    // read white spaces
	    in >> c;
	    while ( c == ' ' ) in >> c;

	    // read '('
	    if ( c != '(' ) {
	      cerr << "QuadraticNumber::operator>> - ( expected" << endl;
	      cout << "got " << c << endl;
	      exit(1);
	    }
	    else  {
	      // read a
	      in >> x.a;

	      // read white spaces
	      in >> c;
	      while ( c == ' ' ) in >> c;

	      // read ','
	      if ( c != ',' ) {
		cerr << "QuadraticNumber::operator>> - , expected" << endl;
		exit(1);
	      }
	      else {
		// read b
		in >> x.b;

		// read white spaces
		in >> c;
		while ( c == ' ' ) in >> c;

		// read ','
		if ( c != ',' ) {
		  cerr << "QuadraticNumber::operator>> - , expected" << endl;
		  exit(1);
		}
		else {
		  // read d
		  in >> x.d;

		  // read white spaces
		  in >> c;
		  while ( c == ' ' ) in >> c;

		  // read ')'
		  if ( c != ')' ) {
		    cerr << "QuadraticNumber::operator>> - ) expected" << endl;
		    exit(1);
		  }
		}
	      }
	    }

	    x.normalize();
	    return in;
      }

      friend std::ostream &
      operator<< (std::ostream & out, const QuadraticNumber<T> & x)
      {
	    out << "(" << x.a << "," << x.b << "," << x.d << ", " << x.QO->getDiscriminant() << ")";
	    return out;
      }
    };



  //
  // Specialized Method Prototypes
  //

  template<>
    void
    QuadraticNumber<GF2EX>::invert ();

  template <> void
  conjugate<GF2EX> (QuadraticNumber<GF2EX> & x,
		    const QuadraticNumber<GF2EX> & y);

  template<>
    void
    mul<GF2EX> (QuadraticNumber<GF2EX> & z, const QuadraticNumber<GF2EX> & x,
		const QuadraticNumber<GF2EX> & y);
  template<>
    void
    sqr<GF2EX> (QuadraticNumber<GF2EX> & z, const QuadraticNumber<GF2EX> & x);

  template<>
    QQ<GF2EX> &
    QuadraticNumber<GF2EX>::getNorm () const;

  template<>
    QQ<GF2EX> &
    QuadraticNumber<GF2EX>::getTrace () const;

} // ANTL

// Unspecialized template definitions.
//#include "../src/Quadratic/QuadraticNumber_impl.hpp"

#endif // ANTL_QUADRATIC_NUMBER_H
