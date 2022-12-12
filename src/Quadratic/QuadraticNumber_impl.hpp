/**
 * @file QuadraticNumber_impl.hpp
 * @author Michael Jacobson
 * @remarks Generic implementation of QuadraticNumber functions.
 */

#include <ANTL/Quadratic/QuadraticNumber.hpp>
namespace NTL {
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
}

namespace ANTL {



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
        return (x.QO->get_discriminant() == y.QO->get_discriminant() && x.a == y.a && x.b = y.b && x.d = y.d);
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
        z.a = x.a * x.d;

        z.b = x.b * x.d;

        z.d = x.a * x.a;
        temp = x.b * x.b;
        temp = temp * x.QO->get_discriminant();
        z.d = z.d - temp;

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

            newA = x.a * y.d;
            temp = y.a * x.d;
            newA = newA + temp;

            newB = x.b * y.d;
            temp = y.b * x.d;
            newB = newB + temp;

            newD = x.d * y.d;

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

            temp = x.d * n;
            z.a = z.a + temp;
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

            newA = x.a * q.getDenominator();
            temp = q.getNumerator() * x.d;
            newA = newA + temp;

            newB = x.b * q.getDenominator();

            newD = x.d * q.getDenominator();

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

            newA = x.a * y.d;
            temp = y.a * x.d;
            newA = newA - temp;

            newB = x.b * y.d;
            temp = y.b * x.d;
            newB = newB - temp;

            newD = x.d * y.d;

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

            temp = x.d * n;
            z.a = z.a - temp;
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

            newA = x.a * q.getDenominator();
            temp = q.getNumerator() * x.d;
            newA = newA - temp;

            newB = x.b * q.getDenominator();

            newD = x.d * q.getDenominator();

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
        if (*z.get_QO() != *x.get_QO() || *z.get_QO() != *y.get_QO()) {
          // TODO:  THROW AN EXCEPTION
        }

        // z = (x.a y.a + x.b y.b N(rho)) + rho(x.a y.b + y.a x.b + x.b y.b Tr(rho))
        //     ---------------------------------------------------------------------
        //                              x.d y.d

        newA = x.get_a() * y.get_a();
        temp = x.get_b() * y.get_b();
        temp = temp * z.get_QO()->get_discriminant();
        newA = newA + temp;

        newB = x.get_a() * y.get_b();
        temp = x.get_b() * y.get_a();
        newB = newB + temp;

        newD = x.get_d() * y.get_d();

        z.set_a(newA);
        z.set_b(newB);
        z.set_d(newD);
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

        z.a = x.a * n;
        z.b = x.b * n;
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

        z.a = x.a * q.getNumerator();
        z.b = x.b * q.getNumerator();
        z.d = z.d * q.getDenominator();
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

        newA = x.a * y.a;
        temp = x.b * y.b;
        temp = temp * z.QO->get_discriminant();
        newA = newA - temp;
        newA = newA * y.d;

        newB = x.b * y.a;
        temp = x.a * y.b;
        newB = newB - temp;
        newB = newB * y.d;

        newD = y.a * y.a;
        temp = y.b * y.b;
        temp = temp * z.QO->get_discriminant();
        newD = temp * x.d;

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
template <class T> void div (QuadraticNumber<T> & z, const QuadraticNumber<T> & x, const T & n) {
  if (z.QO != x.QO) {
    // TODO:  THROW AN EXCEPTION!!!
  }

  z.a = x.a;
  z.b = x.b;
  z.d =  x.d *  n;

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

        z.a = x.a * q.getDenominator();
        z.b = x.b * q.getDenominator();
        z.d = x.d * q.getNumerator();
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

        newA = x.a * x.a;
        temp = x.b * x.b;
        temp = temp * z.QO->get_discriminant();
        newA = newA + temp;

        newB = x.a * x.b;
        newB = newB + newB;

        newD = x.d * x.d;

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
        zz = x * x;
        z = zz * x;
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
        c =  x +  y;
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const T & n)
      {
        QuadraticNumber<T> c;
        c =  x +  n;
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator + (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        QuadraticNumber<T> c;
        c =  x +  q;
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
      QuadraticNumber<T> c;
      c =  x -  y;
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const T & n)
      {
      QuadraticNumber<T> c;
      c =  x -  n;
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator - (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
      QuadraticNumber<T> c;
      c =  x -  q;
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
      QuadraticNumber<T> c;
      c =  x *  y;
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const T & n)
      {
      QuadraticNumber<T> c;
      c =  x *  n;
      return c;
      }

template <class T>
QuadraticNumber<T>
      operator * (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        QuadraticNumber<T> c;
        c =  x *  q;
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & x, const QuadraticNumber<T> & y)
      {
        QuadraticNumber<T> c;
        c =  x / y;
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & a, const T & n)
      {
        QuadraticNumber<T> c;
        c =  a / n;
        return c;
      }

template <class T>
QuadraticNumber<T>
      operator / (const QuadraticNumber<T> & x, const QQ<T> & q)
      {
        QuadraticNumber<T> c;
        c =  x / q;
        return c;
      }

// template <class T>
// const QuadraticNumber<T> &
//       operator += (const QuadraticNumber<T> & x)
//       {
//         *this =  *this +  x;
//         return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator += (const T & n)
//       {
//       *this =  *this +  n;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator += (const QQ<T> & q)
//       {
//       *this =  *this +  q;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator -= (const QuadraticNumber<T> & x)
//       {
//       *this =  *this -  x;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator -= (const T & n)
//       {
//       *this =  *this -  n;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator -= (const QQ<T> & q)
//       {
//       *this =  *this -  q;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator *= (const QuadraticNumber<T> & x)
//       {
//       *this =  *this *  x;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator *= (const T & n)
//       {
//       *this =  *this *  n;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator *= (const QQ<T> & q)
//       {
//       *this =  *this *  q;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator /= (const QuadraticNumber<T> & x)
//       {
//       *this =  *this / x;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator /= (const T & n)
//       {
//       *this =  *this / n;
//       return *this;
//       }
//
// template <class T>
// const QuadraticNumber<T> &
//       operator /= (const QQ<T> & q)
//       {
//       *this =  *this / q;
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
//         out << "(" << x.a << "," << x.b << "," << x.d << ", " << x.QO->get_discriminant() << ")";
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
