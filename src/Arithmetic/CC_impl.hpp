/**
 * @file CC_impl.hpp
 * @remarks This file is to be included from CC.h only.
 * @version $Header$
 *
 * generic methods for CC
 */

namespace ANTL {

  template < class T > CC < T >::CC ()
  {
    //r = i = 0.0;
    ::clear(r);
    ::clear(i);
  }

  template <class T> CC<T>::CC(const CC<T>& x) : r(x.r), i(x.i)
  {
    // r = x.r;
    // i = x.i;
  }

  template <class T> CC<T>::CC(const T& x)
  {
    //r = x;
    //i = 0.0;
    r = x;
    ::clear(i);
  }

  template <class T> CC<T>::CC(const T& rx, const T& ix) : r(rx), i(ix)
  {
    //r = real;
    //i = imag;
  }

  template < class T > CC < T >::~CC ()
  {
  }


  /* access */
  template <class T> void CC<T>::get(T& rx, T& ix)
  {
    rx = r;
    ix = i;
  }

  template <class T> const T& CC<T>::real() const
  {
    return r;
  }

  template <class T> const T& CC<T>::get_r() const
  {
    return r;
  }

  template <class T> const T& CC<T>::imaginary() const
  {
    return i;
  }

  template <class T> const T& CC<T>::get_i() const
  {
    return i;
  }


  /* assignment */
  template <class T> void CC<T>::assign(const T& newreal, const T& newimag)
  {
    r = newreal;
    i = newimag;
  }

  template <class T> void CC<T>::assign(const T& newreal)
  {
    r = newreal;
    ::clear(i);
  }

  template <class T> void CC<T>::assign(const CC<T>& x)
  {
    r = x.r;
    i = x.i;
  }

  template < class T > void CC < T >::clear (void)
  {
    ::clear(r);
    ::clear(i);
  }

  template <class T> void CC<T>::assign_zero()
  {
    ::clear(r);
    ::clear(i);
  }

  template < class T > void CC < T >::set (void)
  {
    ::set(r);
    ::clear(i);
  }

  template <class T> void CC<T>::assign_one()
  {
    ::set(r);
    ::clear(i);
  }

  template < class T > void CC < T >::set_i (void)
  {
    ::clear(r);
    ::set(i);
  }

  template <class T> void swap(CC<T>& x, CC<T>& y)
  {
    CC<T> temp = x;

    x = y;
    y = temp;
  }

  template <class T> const CC<T>& CC<T>::operator = (const CC<T>& A)
  {
    r = A.r;
    i = A.i;
    return *this;
  }

  template <class T> const CC<T>& CC<T>::operator = (const T& A)
  {
    r = A;
    ::clear(i);
    return *this;
  }

  /* comparisions */
  template < class T > bool CC < T >::IsZero()
  {
    return (::IsZero(r) && ::IsZero(i));
  }

  template < class T > bool CC < T >::IsOne()
  {
    return (::IsOne(r) && ::IsZero(i));
  }

  template < class T > bool CC < T >::is_zero()
  {
    return (::IsZero(r) && ::IsZero(i));
  }

  template < class T > bool CC < T >::is_one()
  {
    return (::IsOne(r) && ::IsZero(i));
  }

  template < class T > bool CC < T >::is_i()
  {
    return (::IsZero(r) && ::IsOne(i));
  }

  template < class T > bool CC < T >::is_real (void)
  {
    return (::IsZero(i));
  }

  template <class T> bool CC<T>::is_equal(const CC<T>& B)
  {
    return (r == B.r && i == B.i);
  }

  template <class T> bool CC<T>::is_equal(const T& ea)
  {
    return (r == ea && ::IsZero(i));
  }

  template <class T> bool CC<T>::is_equal(const T& ea, const T& eb)
  {
    return (r == ea && i == eb);
  }

  template <class T> bool operator == (const CC<T>& A, const CC<T>& B)
  {
    return (A.r == B.r && A.i == B.i);
  }

  template <class T> bool operator == (const T& N, const CC<T>& B)
  {
    return (B.r == N && ::IsZero(B.i));
  }

  template <class T> bool operator == (const CC<T>& A, const T& N)
  {
    return (A.r == N && ::IsZero(A.i));
  }

  template <class T> bool operator != (const CC<T>& A, const CC<T>& B)
  {
    return (A.r != B.r || A.i != B.i);
  }

  template <class T> bool operator != (const T& N, const CC<T>& B)
  {
    return (B.r != N || !::IsZero(B.i));
  }

  template <class T> bool operator != (const CC<T>& A, const T& N)
  {
    return (A.r != N || !::IsZero(A.i));
  }


  /* arithmetic */
  template < class T > T CC < T >::norm () const
  {
    return r * r + i * i; // product of number and conjugate
  }

  template <class T> T CC<T>::trace() const
  {
    return r+r; // sum of number and conjugate
  }

  template < class T > void CC < T >::negate (CC < T > &B)
  {
    r = -(B.r);
    i = -(B.i);
  }

  template <class T> void CC<T>::negate()
  {
    r = -r;
    i = -i;
  }

  template < class T > void CC < T >::invert (CC < T > &B)
  {
    T N = B.norm ();

    r = B.r / N;
    i = -B.i / N;
  }

  template <class T> void CC<T>::invert()
  {
    T N = this->norm();

    r = r / N;
    i = -i / N;
  }

  template <class T> void inverse(CC<T>& C, const CC<T>& A)
  {
    T n;

    n = A.norm();
    C.r = A.r / n;
    C.i = -A.i / n;
  }

  template <class T> void conjugate(CC<T>& C, const CC<T>& A)
  {
    C.r = A.r;
    C.i = -A.i;
  }

  /* addition */
  template <class T> void add(CC<T>& C, const CC<T>& A, const CC<T>& B)
  {
    C.r = A.r + B.r;
    C.i = A.i + B.i;
  }

  template <class T> void add(CC<T>& C, const CC<T>& A, const T& B)
  {
    C.r = A.r + B;
    C.i = A.i;
  }

  template <class T> void add(CC<T>& C, const T& A, const T& B)
  {
    C.r = A + B;
    ::clear(C.i);
  }

  /* subtraction */
  template <class T> void subtract(CC<T>& C, const CC<T>& A, const CC<T>& B)
  {
    C.r = A.r - B.r;
    C.i = A.i - B.i;
  }

  template <class T> void subtract(CC<T>& C, const CC<T>& A, const T& B)
  {
    C.r = A.r - B;
    C.i = A.i;
  }

  template <class T> void subtract(CC<T>& C, const T& A, const T& B)
  {
    C.r = A - B;
    ::clear(C.i);
  }

  /* multiplication */
  template <class T> void multiply(CC<T>& C, const CC<T>& A, const CC<T>& B)
  {
    T T1 = (A.r + A.i) * (B.r + B.i);
    T T2 = A.r * B.r;
    T T3 = A.i * B.i;

    C.r = T2 - T3;
    C.i = T1 - T2 - T3;
  }

  template <class T> void multiply(CC<T>& C, const CC<T>& A, const T& B)
  {
    C.r = A.r * B;
    C.i = A.i * B;
  }

  template <class T> void multiply(CC<T>& C, const T& A, const CC<T>& B)
  {
    C.r = B.r * A;
    C.i = B.i * A;
  }

  template <class T> void multiply(CC<T>& C, const T& A, const T& B)
  {
    C.r = A * B;
    ::clear(C.i);
  }

  /* division */
  template <class T> void divide(CC<T>& C, const CC<T>& A, const CC<T>& B)
  {
    T N = B.norm ();

    if (!::IsZero(N))
      {
        T T1 = (A.r + A.i) * (B.r - B.i);
        T T2 = A.r * B.r;
        T T3 = -A.i * B.i;

        C.r = (T2 - T3) / N;
        C.i = (T1 - T2 - T3) / N;
      }
  }

  template <class T> void divide(CC<T>& C, const CC<T>& A, const T& B)
  {
    if (!::IsZero(B))
      {
        C.r = A.r / B;
        C.i = A.i / B;
      }
  }

  template <class T> void divide(CC<T>& C, const T& A, const CC<T>& B)
  {
    T N = B.norm();

    if (!::IsZero(N))
      {
        C.r = A * B.r / N;
        C.i = A * B.i / N;
      }
  }

  template <class T> void divide(CC<T>& C, const T& A, const T& B)
  {
    if (!::IsZero(B))
      {
        C.r = A / B;
        ::clear(C.i);
      }
  }

  /* square */
  template <class T> void square(CC<T>& C, const CC<T>& A)
  {
    multiply(C, A, A);
  }

  template <class T> void square(CC<T>& C, const T& A)
  {
    multiply(C, A, A);
  }

  /* power */
  template <class T> void power(CC<T>& C, const CC<T>& A, const unsigned long e)
  {
    CC<T> temp;

    if (e == 0)
      {
        ::set(C.r);
        ::clear(C.i);
      }
    else if (e == 1)
      {
        C.r = A.r;
        C.i = A.i;
      }
    else if (e > 1)
      {
        if (e & 0x01) // e is odd
	  {
            power(temp, A, (e-1)/2);
            multiply(C, temp, temp);
            multiply(C, C, A);
	  }
        else // e is even
	  {
            power(temp, A, e/2);
            multiply(C, temp, temp);
	  }
      }
  }

  template <class T> void power(CC<T>& C, const CC<T>& A, const ZZ& e)
  {
    ZZ i;
    CC<T> temp;

    if (e == 0)
      {
        ::set(C.r);
        ::clear(C.i);
      }
    else if (e == 1)
      {
        C.r = A.r;
        C.i = A.i;
      }
    else if (e > 1)
      {
        if (bit(e,0)) // e is odd
	  {
            power(temp, A, (e-1)/2);
            multiply(C, temp, temp);
            multiply(C, C, A);
	  }
        else // e is even
	  {
            power(temp, A, e/2);
            multiply(C, temp, temp);
	  }
      }

    /*    if (e > 0)
	  {
	  C = A;
	  i = 1;
	  while (i++ < e)
	  {
	  multiply(C, C, A);
	  }
	  }
    */
  }

  /* overloaded operators */
  template <class T> CC<T> operator - (const CC<T>& A)
  {
    CC<T> C = A;
    C.negate();
    return C;
  }

  template <class T> CC<T> operator + (const CC<T>& A, const CC<T>& B)
  {
    CC<T> C;
    add(C, A, B);
    return C;
  }

  template <class T> CC<T> operator + (const CC<T>& A, const T& B)
  {
    CC<T> C;
    add(C, A, B);
    return C;
  }

  template <class T> CC<T> operator - (const CC<T>& A, const CC<T>& B)
  {
    CC<T> C;
    subtract(C, A, B);
    return C;
  }

  template <class T> CC<T> operator - (const CC<T>& A, const T& B)
  {
    CC<T> C;
    subtract(C, A, B);
    return C;
  }

  template <class T> CC<T> operator * (const CC<T>& A, const CC<T>& B)
  {
    CC<T> C;
    multiply(C, A, B);
    return C;
  }

  template <class T> CC<T> operator * (const CC<T>& A, const T& B)
  {
    CC<T> C;
    multiply(C, A, B);
    return C;
  }

  template <class T> CC<T> operator * (const T& A, const CC<T>& B)
  {
    CC<T> C;
    multiply(C, A, B);
    return C;
  }

  template <class T> CC<T> operator / (const CC<T>& A, const CC<T>& B)
  {
    CC<T> C;
    divide(C, A, B);
    return C;
  }

  template <class T> CC<T> operator / (const CC<T>& A, const T& B)
  {
    CC<T> C;
    divide(C, A, B);
    return C;
  }

  template <class T> CC<T> operator / (const T& A, const CC<T>& B)
  {
    CC<T> C;
    divide(C, A, B);
    return C;
  }

  template <class T> const CC<T>& operator += (CC<T>& A, const CC<T>& B)
  {
    add(A, A, B);
    return A;
  }

  template <class T> const CC<T>& operator -= (CC<T>& A, const CC<T>& B)
  {
    subtract(A, A, B);
    return A;
  }

  template <class T> const CC<T>& operator *= (CC<T>& A, const CC<T>& B)
  {
    multiply(A, A, B);
    return A;
  }

  template <class T> const CC<T>& operator /= (CC<T>& A, const CC<T>& B)
  {
    divide(A, A, B);
    return A;
  }


  /* Input/Output */
  template <class T> std::istream& operator >> (std::istream& in, CC<T>& X)
  {
    long n = 0;
    char c;

    in >> c;
    if (c != '(')
      ANTL_FATAL("ERROR:  CC::operator>>::( expected");

    in >> c;
    while (c != ')' && n != 2)
      {
        in.putback (c);
        if (n == 0)
	  in >> X.r;
        else
	  in >> X.i;
        n++;
        in >> c;
        if (c == ',')
	  in >> c;
      }

    return in;
  }

  template <class T> std::ostream& operator << (std::ostream& out, const CC<T>& X)
  {
    out << "(" << X.r << ", " << X.i << ")" << std::flush;
    return out;
  }


  /* precision */
  template < class T > void CC < T >::SetPrecision (long p)
  {
  }

  template < class T > long CC < T >::Precision (void)
  {
    return 0;
  }

  template < class T > void CC < T >::SetOutputPrecision (long p)
  {
  }

  template < class T > long CC < T >::OutputPrecision (void)
  {
    return 0;
  }

} // ANTL
