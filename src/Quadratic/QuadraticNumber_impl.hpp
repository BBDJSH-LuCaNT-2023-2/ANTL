/**
 * @file QuadraticNumber_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of QuadraticNumber functions.
 */

#include <ANTL/Quadratic/QuadraticNumber.hpp>

namespace ANTL {

template < class T >
void
      clear (QuadraticNumber<T> & z)
      {
        ::clear (z.a);
        ::clear (z.b);
        ::set (z.d);
      }

template < class T >
void
      set (QuadraticNumber<T> & z)
      {
        ::set (z.a);
        ::clear (z.b);
        ::set (z.d);
      }

template < class T >
void
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


template < class T >
bool
      IsZero (const QuadraticNumber<T> & x)
      {
        return (::IsZero (x.a) && ::IsZero (x.b) && ::IsOne (x.d));
      }

template < class T >
bool
      IsOne (const QuadraticNumber<T> & x)
      {
        return (::IsOne (x.a) && ::IsZero (x.b) && ::IsOne (x.d));
      }


template < class T >
bool
      IsEqual (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        return (x.QO == y.QO && x.a == y.a && x.b = y.b && x.d = y.d);
      }

template < class T >
bool
      IsEqual (const QuadraticNumber<T> & x, const T & n)
      {
        return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

template < class T >
bool
      IsEqual (const T & n, const QuadraticNumber<T> & x)
      {
        return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

template < class T >
bool
      IsEqual (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

template < class T >
bool
      IsEqual (const QQ<T> & q, const QuadraticNumber<T> & x)
      {
        return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

template < class T >
bool
      operator == (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        return (x.QO->getDiscriminant() == y.QO->getDiscriminant() && x.a == y.a && x.b = y.b && x.d = y.d);
      }

template < class T >
bool
      operator == (const QuadraticNumber<T> & x, const T & n)
      {
        return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

template < class T >
bool
      operator == (const T & n, const QuadraticNumber<T> & x)
      {
        return (x.a == n && IsZero (x.b) && IsOne (x.d));
      }

template < class T >
bool
      operator == (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

template < class T >
bool
      operator == (const QQ<T> & q, const QuadraticNumber<T> & x)
      {
        return (x.a == q.getNumerator () && IsZero (x.b) && x.d == q.getDenominator ());
      }

template < class T >
bool
      operator != (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        return !(x == y);
      }

template < class T >
bool
      operator != (const QuadraticNumber<T> & x, const T & n)
      {
        return !(x == n);
      }

template < class T >
bool
      operator != (const T & n, const QuadraticNumber<T> & x)
      {
        return !(x == n);
      }

template < class T >
bool
      operator != (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        return !(x == q);
      }

template < class T >
bool
      operator != (const QQ<T> & q, const QuadraticNumber<T> & x)
      {
        return !(x == q);
      }



template < class T>
void
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
template < class T >
void
      inv (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
        static T temp;
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
template <class T>
void
      add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        static T newA, newB, newD, temp;

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
template <class T>
void
      add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n)
      {
        static T temp;

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
template <class T>
void
      add (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        static T newA, newB, newD, temp;
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
template <class T>
void
      sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        static T newA, newB, newD, temp;
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
template < class T>
void
      sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n)
      {
        static T temp;
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
template < class T>
void
      sub (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        static T newA, newB, newD, temp;
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
template <class T>
void
      mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        T newA, newB, newD, temp;
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
template <class T>
void
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
template <class T>
void
      mul (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q)
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
template <class T>
void
      div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        static T newA, newB, newD, temp;
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
template <class T>
void
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
template <class T>
void
      div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const QQ<T> & q)
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
template <class T>
void
      sqr (QuadraticNumber<T> & z, const QuadraticNumber<T> & x)
      {
        static T newA, newB, newD, temp;
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
template <class T>
void
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


template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & a)
      {
        QuadraticNumber<T> c (a);
        c.negate ();
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        QuadraticNumber<T> c;
        add (c, x, y);
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const T & n)
      {
        QuadraticNumber<T> c;
        add (c, x, n);
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        QuadraticNumber<T> c;
        add (c, x, q);
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
      QuadraticNumber<T> c;
      sub (c, x, y);
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const T & n)
      {
      QuadraticNumber<T> c;
      sub (c, x, n);
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
      QuadraticNumber<T> c;
      sub (c, x, q);
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
      QuadraticNumber<T> c;
      mul (c, x, y);
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const T & n)
      {
      QuadraticNumber<T> c;
      mul (c, x, n);
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        QuadraticNumber<T> c;
        mul (c, x, q);
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        QuadraticNumber<T> c;
        div (c, x, y);
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & a, const T & n)
      {
        QuadraticNumber<T> c;
        div (c, a, n);
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        QuadraticNumber<T> c;
        div (c, x, q);
        return c;
      }

// template <class T>
// const QuadraticNumber<T> &
//       operator += (const QuadraticNumber<T> & x)
//       {
//         add (*this, *this, x);
//         return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator += (const T & n)
//       {
//       add (*this, *this, n);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator += (const QQ<T> & q)
//       {
//       add (*this, *this, q);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator -= (const QuadraticNumber<T> & x)
//       {
//       sub (*this, *this, x);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator -= (const T & n)
//       {
//       sub (*this, *this, n);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator -= (const QQ<T> & q)
//       {
//       sub (*this, *this, q);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator *= (const QuadraticNumber<T> & x)
//       {
//       mul (*this, *this, x);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator *= (const T & n)
//       {
//       mul (*this, *this, n);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator *= (const QQ<T> & q)
//       {
//       mul (*this, *this, q);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator /= (const QuadraticNumber<T> & x)
//       {
//       div (*this, *this, x);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator /= (const T & n)
//       {
//       div (*this, *this, n);
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator /= (const QQ<T> & q)
//       {
//       div (*this, *this, q);
//       return *this;
//       }



/**
  *  input / output
  **/

// template <class T>
// std::istream &
//       operator >> <T> (std::istream & in, QuadraticNumber<T> & x)
//       {
//         char c;
//
//         // read white spaces
//         in >> c;
//         while ( c == ' ' ) in >> c;
//
//         // read '('
//         if ( c != '(' ) {
//           cerr << "QuadraticNumber::operator>> - ( expected" << endl;
//           cout << "got " << c << endl;
//           exit(1);
//         }
//         else  {
//           // read a
//           in >> x.a;
//
//           // read white spaces
//           in >> c;
//           while ( c == ' ' ) in >> c;
//
//           // read ','
//           if ( c != ',' ) {
//         cerr << "QuadraticNumber::operator>> - , expected" << endl;
//         exit(1);
//           }
//           else {
//         // read b
//         in >> x.b;
//
//         // read white spaces
//         in >> c;
//         while ( c == ' ' ) in >> c;
//
//         // read ','
//         if ( c != ',' ) {
//           cerr << "QuadraticNumber::operator>> - , expected" << endl;
//           exit(1);
//         }
//         else {
//           // read d
//           in >> x.d;
//
//           // read white spaces
//           in >> c;
//           while ( c == ' ' ) in >> c;
//
//           // read ')'
//           if ( c != ')' ) {
//             cerr << "QuadraticNumber::operator>> - ) expected" << endl;
//             exit(1);
//           }
//         }
//           }
//         }
//
//         x.normalize();
//         return in;
//       }
//
// template <class T>
// std::ostream &
//       operator<<<T> (std::ostream & out, const QuadraticNumber<T> & x)
//       {
//         out << "(" << x.a << "," << x.b << "," << x.d << ", " << x.QO->getDiscriminant() << ")";
//         return out;
//       }

//
// Specialized Method Prototypes
//
template<> void
QuadraticNumber<GF2EX>::invert ();

template <> void
conjugate<GF2EX> (QuadraticNumber<GF2EX> & x, const QuadraticNumber<GF2EX> & y);

template<> void
mul<GF2EX> (QuadraticNumber<GF2EX> & z, const QuadraticNumber<GF2EX> & x, const QuadraticNumber<GF2EX> & y);

template<> void
sqr<GF2EX> (QuadraticNumber<GF2EX> & z, const QuadraticNumber<GF2EX> & x);

// template<> QQ<GF2EX> &
// QuadraticNumber<GF2EX>::getNorm () const;

template<> QQ<GF2EX> &
QuadraticNumber<GF2EX>::getTrace () const;

template<class T> bool
QuadraticNumber<T>::isUnit () const
{
  // Note that in the function field case, the normalization we have imposed implies that if this is
  // a unit, then the norm will be equal to 1 (not just a constant)

  return deg(getNorm()) == 0;
}

} // ANTL
