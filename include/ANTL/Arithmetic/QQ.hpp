/**
 * @file QQ.hpp
 * @author Jonathan Hammell / Michael Jacobson
 * @brief template rational number class.  The numerator and denominator are always
 * stored in lowest terms.
 */

#ifndef ANTL_QQ_H
#define ANTL_QQ_H

#include <ANTL/common.hpp>

namespace ANTL
{
  // class prototype needed by friend function definitions
  template<class T>
    class QQ;

  // Forward Declarations of friend functions
  template<class T>
    void
    clear (QQ<T> & x);

  template<class T>
    void
    set (QQ<T> & x);

  template<class T>
    void
    assign (QQ<T> & z, const QQ<T> & x);

  template<class T>
    bool
    IsZero (const QQ<T>& x);

  template<class T>
    bool
    IsOne (const QQ<T>& x);

  template<class T>
    bool
    IsEqual (const QQ<T> & x, const QQ<T> & y);

  template<class T>
    bool
    IsEqual (const QQ<T> & x, const T & n);

  template<class T>
    bool
    IsEqual (const T & n, const QQ<T> & x);



  template<class T>
    bool
    operator == (const QQ<T> & x, const QQ<T> & y);

  template<class T>
    bool
    operator == (const QQ<T> & x, const T & n);

  template<class T>
    bool
    operator == (const T & n, const QQ<T> & x);

  template<class T>
    bool
    operator != (const QQ<T> & x, const QQ<T> & y);

  template<class T>
    bool
    operator != (const QQ<T> & x, const T & n);

  template<class T>
    bool
    operator != (const T & n, const QQ<T> & x);

  template<class T>
    void
    inv (QQ<T> & x, const QQ<T> & y);

  template<class T>
    void
    add (QQ<T> & z, const QQ<T> & x, const QQ<T> & y);

  template<class T>
    void
    add (QQ<T> & z, const QQ<T> & x, const T & n);

  template<class T>
    void
    add (QQ<T> & z, const T & n, const QQ<T> & x);

  template<class T>
    void
    sub (QQ<T> & z, const QQ<T> & x, const QQ<T> & y);

  template<class T>
    void
    sub (QQ<T> & z, const QQ<T> & x, const T & n);

  template<class T>
    void
    sub (QQ<T> & z, const T & n, const QQ<T> & x);

  template<class T>
    void
    mul (QQ<T> & z, const QQ<T> & x, const QQ<T> & y);

  template<class T>
    void
    mul (QQ<T> & z, const QQ<T> & x, const T & n);

  template<class T>
    void
    mul (QQ<T> & z, const T & n, const QQ<T> & x);

  template<class T>
    void
    div (QQ<T> & z, const QQ<T> & x, const QQ<T> & y);

  template<class T>
    void
    div (QQ<T> & z, const QQ<T> & x, const T & n);

  template<class T>
    void
    div (QQ<T> & z, const T & n, const QQ<T> & x);

  template<class T>
    void
    sqr (QQ<T> & z, const QQ<T> & x);

  template<class T>
    void
    cube (QQ<T> & z, const QQ<T> & x);

  template<class T>
    QQ<T>
    operator - (const QQ<T> & x);

  template<class T>
    QQ<T>
    operator + (const QQ<T> & x, const QQ<T> & y);

  template<class T>
    QQ<T>
    operator + (const QQ<T> & x, const T & n);

  template<class T>
    QQ<T>
    operator + (const T & n, const QQ<T> & x);

  template<class T>
    QQ<T>
    operator - (const QQ<T> & x, const QQ<T> & b);

  template<class T>
    QQ<T>
    operator - (const QQ<T> & x, const T & n);

  template<class T>
    QQ<T>
    operator - (const T &, const QQ<T> & x);

  template<class T>
    QQ<T>
    operator * (const QQ<T> & x, const QQ<T> & y);

  template<class T>
    QQ<T>
    operator * (const QQ<T> & x, const T & n);

  template<class T>
    QQ<T>
    operator * (const T & n, const QQ<T> & x);

  template<class T>
    QQ<T>
    operator / (const QQ<T> & x, const QQ<T> & y);

  template<class T>
    QQ<T>
    operator / (const QQ<T> & x, const T & n);

  template<class T>
    QQ<T>
    operator / (const T & n, const QQ<T> & x);

  template<class T>
    std::istream &
    operator>> (std::istream & in, QQ<T> & x);

  template<class T>
    std::ostream &
    operator<< (std::ostream & out, const QQ<T> & x);

  /**
   * @brief Rational number or function
   * @remarks An element of QQ represents either a rational number or a rational function over a finite field.
   *
   * Instances of QQ are all normalized such that the numerator and denominator are relatively prime.  In the
   * case of rational functions, the denominator is monic.
   *
   * Currently, the following template parameters T are instantiated:
   *    base types:
   *       long --- rational number with word-sized coefficients
   *    NTL:
   *       ZZ --- rational number (arbitrary sized coefficients)
   *       ZZ_pX --- rational function over Fp (char <> 2)
   *       ZZ_pEX --- rational function over Fq (char <> 2)
   *       zz_pX --- rational function over Fp (char <> 2, p < 2^64)
   *       zz_pEX --- rational function over Fq (char <> 2, p < 2^64)
   *       GF2EX --- rational function over Fq (char = 2)
   */
  template<class T>
    class QQ
    {
    protected:
      T a; /**< numerator */
      T d; /**< denominator */

      static T newA, newD, temp;// temporary variables for arithmetic operations

      void
      normalize ()
      {
	// remove common factors
	T g = GCD (a, d);
	if (!::IsOne (g))
	  {
	    ::divide (a, a, g);
	    ::divide (d, d, g);
	  }
        //make monic
        if (!::IsOne(LeadCoeff(d)))
        {
          NTL::div(a, a, LeadCoeff(d));
          NTL::div(d, d, LeadCoeff(d));
        }
      }

    public:
      QQ ()
      {
	::clear (a);
	NTL::set (d);
      }

      /**
       * @brief Copy constructor
       * @param[in] x the QQ to copy
       */
      QQ (const QQ<T>& x)
      {
	a = x.a;
	d = x.d;
      }

      /**
       * @brief Copy constructor
       * @param[in] x the constant to copy
       */
      QQ (const T& x)
      {
	a = x;
	NTL::set (d);
      }

      /**
       * @brief Constructs rational number from the given coefficients
       * @param[in] num the numerator
       * @param[in] den the denominator
       */
      QQ (const T& num, const T& den)
      {
	a = num;
	d = den;
	normalize ();
      }

      ~QQ ()
      {

      }

      /**
       * Accessor methods
       */

      const T&
      getNumerator () const
      {
	return a;
      }

      const T&
      getN () const
      {
	return a;
      }

      const T&
      getDenominator () const
      {
	return d;
      }

      const T&
      getD () const
      {
	return d;
      }

      /**
       *  Assignments and mutator methods.
       **/

      void
      setNumerator (const T & na)
      {
	a = na;
	normalize ();
      }

      void
      setN (const T & na)
      {
	a = na;
	normalize ();
      }

      void
      setDenominator (const T & den)
      {
	d = den;
	normalize ();
      }

      void
      setD (const T & den)
      {
	d = den;
	normalize ();
      }

      /**
       * @brief Sets z to zero
       * @param[out] z QQ to be set to zero
       */
      friend void
      clear (QQ<T> & z)
      {
	::clear (z.a);
	NTL::set (z.d);
      }

      /**
       * @brief Sets z to one
       * @param[out] z QQ to be set to one
       */
      friend void
      set (QQ<T> & z)
      {
	NTL::set (z.a);
	NTL::set (z.d);
      }

      /**
       * @brief Sets this QQ equal to x
       * @param[in] x value to give the QQ
       */
      void
      assign (const QQ<T> & x)
      {
	a = x.a;
	d = x.d;
      }

      /**
       * @brief Sets this QQ equal to n
       * @param[in] n value to give the QQ
       */
      void
      assign (const T & n)
      {
	a = n;
	NTL::set (d);
      }

      /**
       * @brief Sets the QQ z equal to x
       * @param[out] z result of assignment
       * @param[in] x value to be assigned
       */
      friend void
      assign (QQ<T> & z, const QQ<T> & x)
      {
	z.a = x.a;
	z.d = x.d;
      }

      /**
       * @brief Assignment operator
       * @param[in] x value to give the QQ
       */
      const QQ<T> &
      operator= (const QQ<T> & x)
      {
	this->assign (x);
	return *this;

      }

      /**
       * Comparison methods
       **/

      /**
       * @brief Tests whether this QQ is equal to zero
       * @return True if the QQ is zero
       */
      bool
      isZero () const
      {
	return (::IsZero (a) && ::IsOne (d));
      }
      ;

      /**
       * @brief Tests whether this QQ is equal to one
       * @return True if the QQ is one
       */
      bool
      isOne () const
      {
	return (::IsOne (a) && ::IsOne (d));
      }

      /**
       * @brief Tests whether x is equal to zero
       * @param[in] x QQ to test
       * @return True if x is zero
       */
      friend bool
      IsZero (const QQ<T> & x)
      {
	return (::IsZero (x.a) && ::IsOne (x.d));
      }

      /**
       * @brief Tests whether x is equal to one
       * @param[in] x QQ to test
       * @return True if x is one
       */
      friend bool
      IsOne (const QQ<T> & x)
      {
	return (::IsOne (x.a) && ::IsOne (x.d));
      }

      /**
       * @brief Determines whether this QQ is an integer (or polynomial)
       * @return true if this is an integer (or polynomial)
       */
      bool
      isInteger () const
      {
	return ::IsOne (d);
      }

      /**
       * @brief Tests whether this QQ is equal to y
       * @param[in] y QQ to compare
       * @return True if the QQ is equal to y
       */
      bool
      isEqual (const QQ<T> & x)
      {
	return (a == x.a && d == x.d);
      }

      /**
       * @brief Tests whether this QQ is equal to n
       * @param[in] n constant to compare
       * @return True if the QQ is equal to n
       */
      bool
      isEqual (const T & n)
      {
	return (a == n && ::IsOne (d));
      }

      /**
       * @brief Tests whether x is equal to y
       * @return True if x is equal to y
       */
      friend bool
      IsEqual (const QQ<T> & x, const QQ<T> & y)
      {
	return (x.a == y.a && x.d == y.d);
      }

      /**
       * @brief Tests whether x is equal to the constant n
       * @param[in] x QQ to compare
       * @param[in] n constant to compare
       * @return True if x is equal to n
       */
      friend bool
      IsEqual (const QQ<T> & x, const T & n)
      {
	return (x.a == n && ::IsOne (x.d));
      }

      /**
       * @brief Tests whether x is equal to the constant n
       * @param[in] n constant to compare
       * @param[in] x QQ to compare
       * @return True if x is equal to n
       */
      friend bool
      IsEqual (const T & n, const QQ<T> & x)
      {
	return (x.a == n && ::IsOne (x.d));
      }

      /**
       * @brief Equality test operator
       * @param[in] x QQ to compare
       * @param[in] y QQ to compare
       * @return True if x is equal to y
       */
      friend bool
      operator == (const QQ<T> & x, const QQ<T> & y)
      {
	return (x.a == y.a && x.d == y.d);
      }

      /**
       * @brief Equality test operator
       * @param[in] x QQ to compare
       * @param[in] n constant to compare
       * @return True if x is equal to n
       */
      friend bool
      operator == (const QQ<T> & x, const T & n)
      {
	return (x.a == n && ::IsOne (x.d));
      }

      /**
       * @brief Equality test operator
       * @param[in] n constant to compare
       * @param[in] x QQ to compare
       * @return True if x is equal to n
       */
      friend bool
      operator == (const T & n, const QQ<T> & x)
      {
	return (x.a == n && ::IsOne (x.d));
      }

      /**
       * @brief Inequality test operator
       * @param[in] x QQ to compare
       * @param[in] y QQ to compare
       * @return True if x is not equal to y
       */

      friend bool
      operator != (const QQ<T> & x, const QQ<T> & y)
      {
	return (x.a != y.a || x.d != y.d);
      }

      /**
       * @brief Inequality test operator
       * @param[in] x QQ to compare
       * @param[in] n constant to compare
       * @return True if x is not equal to n
       */
      friend bool
      operator != (const QQ<T> & x, const T & n)
      {
	return (x.a != n || !::IsOne (x.d));
      }

      /**
       * @brief Inequality test operator
       * @param[in] n constant to compare
       * @param[in] x QQ to compare
       * @return True if x is not equal to n
       */
      friend bool
      operator != (const T & n, const QQ<T> & x)
      {
	return (x.a != n || !::IsOne (x.d));
      }

      /**
       * Arithmetic operations
       **/

      /**
       * @brief Negates the QQ
       */
      void
      negate ()
      {
	a = -a;
      }

      /**
       * @brief Inverts the QQ
       */
      void
      invert ()
      {

	temp = d;
	d = a;
	a = temp;
        normalize();
      }

      /**
       * @brief Sets the QQ z equal to the inverse of x
       * @param[out] z inverse of x
       * @param[in] x value to be inverted
       */
      friend void
      inv (QQ<T> & z, const QQ<T> & x)
      {
	z.a = x.d;
	z.d = x.a;
        z.normalize();
      }

      /**
       * @brief Computes the sum of x and y
       * @param[out] z = x + y
       * @param[in] x first summand
       * @param[in] y second summand
       */
      friend void
      add (QQ<T> & z, const QQ<T> & x, const QQ<T> & y)
      {
	if (x.isZero ())
	  z.assign (y);
	else if (y.isZero ())
	  z.assign (x);
	else
	  {
	    //     y.d * x.a + x.d * y.a
	    // z = ---------------------
	    //          x.d * y.d

	    ::mul (newA, x.a, y.d);
	    ::mul (temp, y.a, x.d);
	    ::add (newA, newA, temp);

	    ::mul (newD, x.d, y.d);

	    z.a = newA;
	    z.d = newD;
	    z.normalize ();
	  }
      }


      /**
       * @brief Computes the sum of x and n (a constant)
       * @param[out] z = x + n
       * @param[in] x first summand
       * @param[in] n second summand (a constant)
       */
      friend void
      add (QQ<T> & z, const QQ<T> & x, const T & n)
      {
	if (x.isZero ())
	  z.assign (n);
	else if (::IsZero (n))
	  z.assign (x);
	else
	  {
	    //     x.a + x.d * n
	    // z = -------------
	    //         x.d

	    ::mul (temp, x.d, n);
	    ::add (z.a, x.a, temp);
	    z.d = x.d;
	    z.normalize ();
	  }

      }


      /**
       * @brief Computes the sum of x and n (a constant)
       * @param[out] z = x + n
       * @param[in] n first summand (a constant)
       * @param[in] x second summand
       */
      friend void
      add (QQ<T> & z, const T & n, const QQ<T> & x)
      {
	if (x.isZero ())
	  z.assign (n);
	else if (::IsZero (n))
	  z.assign (x);
	else
	  {
	    //     x.a + x.d * n
	    // z = -------------
	    //         x.d

	    ::mul (temp, x.d, n);
	    ::add (z.a, x.a, temp);
	    z.d = x.d;
	    z.normalize ();
	  }

      }

      /**
       * @brief Computes the difference of x and y
       * @param[out] z = x - y
       * @param[in] x first term
       * @param[in] y second term
       */
      friend void
      sub (QQ<T> & z, const QQ<T> & x, const QQ<T> & y)
      {
	if (x.isZero ())
	  {
	    z.assign (y);
	    z.negate ();
	  }
	else if (y.isZero ())
	  z.assign (x);
	else
	  {
	    //     y.d * x.a - x.d * y.a
	    // z = ---------------------
	    //          x.d * y.d

	    ::mul (newA, x.a, y.d);
	    ::mul (temp, y.a, x.d);
	    ::sub (newA, newA, temp);

	    ::mul (newD, x.d, y.d);

	    z.a = newA;
	    z.d = newD;
	    z.normalize ();
	  }
      }

      /**
       * @brief Computes the difference of x and n (a constant)
       * @param[out] z = x - n
       * @param[in] x first term
       * @param[in] n second term (a constant)
       */
      friend void
      sub (QQ<T> & z, const QQ<T> & x, const T & n)
      {
	if (x.isZero ())
	  {
	    z.assign (n);
	    z.negate ();
	  }
	else if (::IsZero (n))
	  z.assign (x);
	else
	  {
	    //     x.a - x.d * n
	    // z = -------------
	    //         x.d

	    ::mul (temp, x.d, n);
	    ::sub (z.a, x.a, temp);
	    z.d = x.d;
	    z.normalize ();
	  }
      }

      /**
       * @brief Computes the difference of x and n (a constant)
       * @param[out] z = n - x
       * @param[in] n first term (a constant)
       * @param[in] x second term
       */
      friend void
      sub (QQ<T> & z, const T & n, const QQ<T> & x)
      {
	if (x.isZero ())
	  {
	    z.assign (n);
	  }
	else if (::IsZero (n))
         {
	  z.assign (x);
          z.negate();
         }
	else
	  {
	    //     x.d * n - x.a
	    // z = -------------
	    //         x.d

	    ::mul (temp, x.d, n);
	    ::sub (z.a, temp, x.a);
	    z.d = x.d;
	    z.normalize ();
	  }
      }

      /**
       * @brief Computes the product of x and y
       * @param[out] z = x * y
       * @param[in] x first term
       * @param[in] y second term
       */
      friend void
      mul (QQ<T> & z, const QQ<T> & x, const QQ<T> & y)
      {
	::mul (z.a, x.a, y.a);
	::mul (z.d, x.d, y.d);
	z.normalize ();
      }

      /**
       * @brief Computes the product of x and n (a constant)
       * @param[out] z = x * n
       * @param[in] x first term
       * @param[in] n second term (a constant)
       */
      friend void
      mul (QQ<T> & z, const QQ<T> & x, const T & n)
      {
	::mul (z.a, x.a, n);
	z.d = x.d;
	z.normalize ();
      }

      /**
       * @brief Computes the product of x and n (a constant)
       * @param[out] z = x * n
       * @param[in] n first term (a constant)
       * @param[in] x second term
       */
      friend void
      mul (QQ<T> & z, const T & n, const QQ<T> & x)
      {
	::mul (z.a, x.a, n);
	z.d = x.d;
	z.normalize ();
      }

      /**
       * @brief Computes the quotient of x and y
       * @param[out] z = x / y
       * @param[in] x first term
       * @param[in] y second term
       */
      friend void
      div (QQ<T> & z, const QQ<T> & x, const QQ<T> & y)
      {
	::mul (z.a, x.a, y.d);
	::mul (z.d, x.d, y.a);
	z.normalize ();
      }

      /**
       * @brief Computes the quotient of x and n (a constant)
       * @param[out] z = x / n
       * @param[in] x first term
       * @param[in] n second term (a constant)
       */
      friend void
      div (QQ<T> & z, const QQ<T> & x, const T & n)
      {
	z.a = x.a;
	::mul (z.d, x.d, n);
	z.normalize ();
      }

      /**
       * @brief Computes the quotient of n (a constant) and x
       * @param[out] z = n / x
       * @param[in] x first term
       * @param[in] n second term (a constant)
       */
      friend void
      div (QQ<T> & z, const T & n, const QQ<T> & x)
      {
	::mul (z.a, x.d, n);
	z.d = x.a;
	z.normalize ();
      }

      /**
       * @brief Computes the square of x
       * @param[out] z = x^2
       * @param[in] x value to square
       */
      friend void
      sqr (QQ<T> & z, const QQ<T> & x)
      {
	::sqr (z.a, x.a);
	::sqr (z.d, x.d);
      }

      /**
       * @brief Computes the cube of x
       * @param[out] z = x^3
       * @param[in] x value to cube
       */
      friend void
      cube (QQ<T> & z, const QQ<T> & x)
      {
	::sqr (newA, x.a);
	::mul (newA, newA, x.a);

	::sqr (newD, x.d);
	::mul (newD, newD, x.d);

	z.a = newA;
	z.d = newD;
      }

      friend QQ<T>
      operator - (const QQ<T> & x)
      {
	QQ<T> c (x);
	c.negate ();
	return c;
      }

      friend QQ<T>
      operator + (const QQ<T> & x, const QQ<T> & y)
      {
	QQ<T> c;
	add (c, x, y);
	return c;
      }

      friend QQ<T>
      operator + (const QQ<T> & x, const T & n)
      {
	QQ<T> c;
	add (c, x, n);
	return c;
      }

      friend QQ<T>
      operator + (const T & n, const QQ<T> & x)
      {
	QQ<T> c;
	add (c, n, x);
	return c;
      }

      friend QQ<T>
      operator - (const QQ<T> & x, const QQ<T> & y)
      {
	QQ<T> c;
	sub (c, x, y);
	return c;
      }

      friend QQ<T>
      operator - (const QQ<T> & x, const T & n)
      {
	QQ<T> c;
	sub (c, x, n);
	return c;
      }

      friend QQ<T>
      operator - (const T & n, const QQ<T> & x)
      {
	QQ<T> c;
	sub (c, n, x);
	return c;
      }

      friend QQ<T>
      operator * (const QQ<T> & x, const QQ<T> & y)
      {
	QQ<T> c;
	mul (c, x, y);
	return c;
      }

      friend QQ<T>
      operator * (const QQ<T> & x, const T & n)
      {
	QQ<T> c;
	mul (c, x, n);
	return c;
      }

      friend QQ<T>
      operator * (const T & n, const QQ<T> & x)
      {
	QQ<T> c;
	mul (c, n, x);
	return c;
      }

      friend QQ<T>
      operator / (const QQ<T> & x, const QQ<T> & y)
      {
	QQ<T> c;
	div (c, x, y);
	return c;
      }

      friend QQ<T>
      operator / (const QQ<T> & x, const T & n)
      {
	QQ<T> c;
	div (c, x, n);
	return c;
      }

      friend QQ<T>
      operator / (const T & n, const QQ<T> & x)
      {
	QQ<T> c;
	div (c, n, x);
	return c;
      }

      const QQ<T> &
      operator += (const QQ<T> & x)
      {
	add (*this, *this, x);
	return *this;
      }

      const QQ<T> &
      operator += (const T & n)
      {
	add (*this, *this, n);
	return *this;
      }

      const QQ<T> &
      operator -= (const QQ<T> & x)
      {
	sub (*this, *this, x);
	return *this;
      }

      const QQ<T> &
      operator -= (const T & n)
      {
	sub (*this, *this, n);
	return *this;
      }

      const QQ<T> &
      operator *= (const QQ<T> & x)
      {
	mul (*this, *this, x);
	return *this;
      }

      const QQ<T> &
      operator *= (const T & n)
      {
	mul (*this, *this, n);
	return *this;
      }

      const QQ<T> &
      operator /= (const QQ<T> & x)
      {
	div (*this, *this, x);
	return *this;
      }

      const QQ<T> &
      operator /= (const T & n)
      {
	div (*this, *this, n);
	return *this;
      }

      /**
       *  input / output
       **/

      friend std::istream &
      operator>> (std::istream & in, QQ<T> & x)
      {
	char c;

	// read white spaces
	in >> c;
	while (c == ' ')
	  in >> c;

	// read '('
	if (c != '(')
	  {
	    cerr << "QQ::operator>> - ( expected" << endl;
	    cout << "got " << c << endl;
	    exit (1);
	  }
	else
	  {
	    // read a
	    in >> x.a;

	    // read white spaces
	    in >> c;
	    while (c == ' ')
	      in >> c;

	    // read ','
	    if (c != ',')
	      {
		cerr << "QQ::operator>> - , expected" << endl;
		exit (1);
	      }
	    else
	      {
		// read d
		in >> x.d;

		// read white spaces
		in >> c;
		while (c == ' ')
		  in >> c;

		// read ')'
		if (c != ')')
		  {
		    cerr << "QQ::operator>> - ) expected" << endl;
		    exit (1);
		  }
	      }
	  }

	x.normalize ();
	return in;
      }

      friend std::ostream &
      operator<< (std::ostream & out, const QQ<T> & x)
      {
	out << "(" << x.a << "," << x.d << ")";
	return out;
      }
    };

  //static member initializations
  template<class T>
  T QQ<T>::temp;

  template<class T>
  T QQ<T>::newA;

  template<class T>
  T QQ<T>::newD;


  // Partial specializations
  template<>
    inline void
    QQ<long>::normalize (){

       long g = GCD(this->a, this->d);
       if(g > 1){
           this->a /= g;
           this->d /= g;
       }
       if(this->d < 0){
          this->a = -this->a;
          this->d = -this->d;
       }
    };

  template<>
    inline void
    QQ<ZZ>::normalize (){

       ZZ g = GCD(this->a, this->d);
       if(g > 1){
           this->a /= g;
           this->d /= g;
       }
       if(this->d < 0){
          this->a = -this->a;
          this->d = -this->d;
       }
    };


} // ANTL

#endif // ANTL_QQ_H
