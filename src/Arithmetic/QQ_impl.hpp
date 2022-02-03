/**
 * @file QQ_impl.hpp
 * @remarks This file is to be included from QQ.hpp only.
 * @version $Header$
 *
 * generic methods for QQ
 */

namespace ANTL {

  template <class T> QQ<T>::QQ() {
    ::clear(n);
    ::set(d);
  }

  template <class T>QQ<T>::QQ(const QQ<T>& x) : n(x.n), d(x.d) {
    //n = x.n;
    //d = x.d;
  }

  template <class T> QQ<T>::QQ(const T& x) {
    n = x;
    ::set(d);
  }

  template <class T> QQ<T>::QQ(const T& num, const T& den) {
    if (::IsZero(den)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    n = num;
    d = den;
    normalize();
  }

  template <class T> QQ<T>::~QQ() {}


  /* access */
  template <class T> void QQ<T>::get(T& num, T& den) {
    num = n;
    den = d;
  }

  template <class T> const T& QQ<T>::numerator() const {
    return n;
  }

  template <class T> const T& QQ<T>::get_n() const {
    return n;
  }

  template <class T> const T& QQ<T>::denominator() const {
    return d;
  }

  template <class T> const T& QQ<T>::get_d() const {
    return d;
  }


  /* assignment */
  template <class T> void QQ<T>::assign(const T& newnum, const T& newden) {
    if (::IsZero(newden)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    n = newnum;
    d = newden;
    normalize();
  }

  template <class T> void QQ<T>::assign(const T& newnum) {
    n = newnum;
    ::set(d);
  }

  template <class T> void QQ<T>::assign(const QQ<T>& x) {
    n = x.n;
    d = x.d;
  }

  template <class T> void QQ<T>::clear(void) {
    ::clear(n);
    ::set(d);
  }

  template <class T> void QQ<T>::assign_zero() {
    ::clear(n);
    ::set(d);
  }

  template <class T> void QQ<T>::set(void) {
    ::set(n);
    ::set(d);
  }

  template <class T> void QQ<T>::assign_one() {
    ::set(n);
    ::set(d);
  }

  template <class T> void swap(QQ<T>& x, QQ<T>& y) {
    QQ<T> temp = x;
    x = y;
    y = temp;
  }

  template <class T> const QQ<T>& QQ<T>::operator = (const QQ<T>& A) {
    n = A.n;
    d = A.d;
    return *this;
  }

  template <class T> const QQ<T>& QQ<T>::operator = (const T& A) {
    n = A;
    ::set(d);
    return *this;
  }

  /* comparisions */
  template <class T> bool QQ<T>::is_zero() {
    return (::IsZero(n));
  }

  template <class T> bool QQ<T>::is_one() {
    return (::IsOne(n) && ::IsOne(d));
  }

  template <class T> bool IsZero(const QQ<T>& A) {
    return (::IsZero(A.n)); 
  }

  template <class T> bool IsOne(const QQ<T>& A) {
    return (::IsOne(A.n) && ::IsOne(A.d)); 
  }

  template < class T > bool QQ < T >::is_integer(void) {
    return (::IsOne(d));
  }

  template <class T> bool QQ<T>::is_equal(const QQ<T>& B) {
    return (n == B.n && d == B.d);
  }

  template <class T> bool QQ<T>::is_equal(const T& ea) {
    return (n == ea && ::IsOne(d));
  }

  template <class T> bool QQ<T>::is_equal(const T& ea, const T& eb) {
    return (n == ea && d == eb);
  }

  template <class T> bool operator == (const QQ<T>& A, const QQ<T>& B) {
    return (A.n == B.n && A.d == B.d);
  }

  template <class T> bool operator == (const T& N, const QQ<T>& B) {
    return (B.n == N && ::IsOne(B.d));
  }

  template <class T> bool operator == (const QQ<T>& A, const T& N) {
    return (A.n == N && ::IsOne(A.d));
  }

  template <class T> bool operator != (const QQ<T>& A, const QQ<T>& B) {
    return (A.n != B.n || A.d != B.d);
  }

  template <class T> bool operator != (const T& N, const QQ<T>& B) {
    // Reduce?
    return (B.n != N || ::IsOne(B.d));
  }

  template <class T> bool operator != (const QQ<T>& A, const T& N) {
    // Reduce?
    return (A.n != N || !::IsOne(A.d));
  }


  /* arithmetic */
  template <class T> void QQ<T>::negate(const QQ<T> &B) {
    n = -(B.n);
  }

  template <class T> void QQ<T>::negate() {
    n = -n;
  }

  template <class T> void QQ<T>::invert(const QQ <T> &B) {
    if (::IsZero(B.n)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    n = B.d;
    d = B.n;
    normalize();
  }

  template <class T> void QQ<T>::invert() {
    if (::IsZero(n)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    T temp = n;
    n = d;
    d = temp;
    normalize();
  }

  /*template <class T> void inverse(QQ<T>& C, const QQ<T>& A) {
    T n;

    n = A.norm();
    C.r = A.r / n;
    C.i = -A.i / n;
    }*/

  /* addition */
  template <class T> void add(QQ<T>& C, const QQ<T>& A, const QQ<T>& B) {
    if (A.d == B.d) {
      C.n = A.n + B.n;
      C.d = A.d;
    } else {
      // We could compute LCM, but probably not needed 
      // unless used for really big numbers
      C.n = (A.n * B.d) + (B.n * A.d);
      C.d = A.d * B.d;
    }
    C.normalize();
  }

  template <class T> void add(QQ<T>& C, const QQ<T>& A, const T& B) {
    T Bn = B * A.d;
    C.n = A.n + Bn;
    C.d = A.d;
    C.normalize();
  }

  template <class T> void add(QQ<T>& C, const T& A, const T& B) {
    C.n = A + B;
    ::set(C.d);
  }

  /* subtraction */
  template <class T> void subtract(QQ<T>& C, const QQ<T>& A, const QQ<T>& B) {
    if (A.d == B.d) {
      C.n = A.n - B.n;
      C.d = A.d;
    } else {
      // We could compute LCM, but probably not needed 
      // unless used for really big numbers
      C.n = (A.n * B.d) - (B.n * A.d);
      C.d = A.d * B.d;
    }
    C.normalize();
  }

  template <class T> void subtract(QQ<T>& C, const QQ<T>& A, const T& B) {
    T Bn = B * A.d;
    C.n = A.n - Bn;
    C.d = A.d;
    C.normalize();
  }

  template <class T> void subtract(QQ<T>& C, const T& A, const T& B) {
    C.n = A - B;
    ::set(C.d);
  }

  /* multiplication */
  template <class T> void multiply(QQ<T>& C, const QQ<T>& A, const QQ<T>& B) {
    // We could possibly reduce earlier, but probably not needed 
    // unless used for really big numbers
    C.n = A.n * B.n;
    C.d = A.d * B.d;
    C.normalize();
  }

  template <class T> void multiply(QQ<T>& C, const QQ<T>& A, const T& B) {
    // We could possibly reduce earlier, but probably not needed 
    // unless used for really big numbers
    C.n = A.n * B;
    C.d = A.d;
    C.normalize();
  }

  template <class T> void multiply(QQ<T>& C, const T& A, const QQ<T>& B) {
    // We could possibly reduce earlier, but probably not needed 
    // unless used for really big numbers
    C.n = A * B.n;
    C.d = B.d;
    C.normalize();
  }

  template <class T> void multiply(QQ<T>& C, const T& A, const T& B) {
    C.n = A * B;
    ::set(C.d);
  }

  /* division */
  template <class T> void divide(QQ<T>& C, const QQ<T>& A, const QQ<T>& B) {
    if (::IsZero(B.n)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    // We could possibly reduce earlier, but probably not needed 
    // unless used for really big numbers
    C.n = A.n * B.d;
    C.d = A.d * B.n;
    C.normalize();
  }

  template <class T> void divide(QQ<T>& C, const QQ<T>& A, const T& B) {
    if (::IsZero(B)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    // We could possibly reduce earlier, but probably not needed 
    // unless used for really big numbers
    C.n = A.n;
    C.d = A.d * B;
    C.normalize();
  }

  template <class T> void divide(QQ<T>& C, const T& A, const QQ<T>& B) {
    if (::IsZero(B.n)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    // We could possibly reduce earlier, but probably not needed 
    // unless used for really big numbers
    C.n = A * B.d;
    C.d = B.n;
    C.normalize();
  }

  template <class T> void divide(QQ<T>& C, const T& A, const T& B) {
    if (::IsZero(B)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
    C.n = A;
    C.d = B;
    C.normalize();
  }

  /* square */
  template <class T> void square(QQ<T>& C, const QQ<T>& A) {
    multiply(C, A, A);
  }

  template <class T> void square(QQ<T>& C, const T& A) {
    multiply(C, A, A);
  }

  /* power */
  template <class T> void power(QQ<T>& C, const QQ<T>& A, unsigned long e) {
    if (e == 0) {
      C.assign_one();
    } else if (e == 1) {
      C.assign(A);
    } else if (e > 1) {
      QQ<T> B;
      B.assign(A);
      C.assign_one();
      while (e > 0) {
	if (e & 1)
	  multiply(C, C, B);
	e >>= 1;
	if (e > 0)
	  square(B, B);
      }
    }
  }

  template <class T> void power(QQ<T>& C, const QQ<T>& A, long e) {
    unsigned long exp;
    exp = abs(e);
	
    if (e < 0) {
      QQ<T> B;
      B.invert(A);
      power(C, B, exp);
    } else {
      power(C, A, exp);
    }
  }

  template <class T> void power(QQ<T>& C, const QQ<T>& A, const ZZ& e) {
    if (e == 0) {
      C.assign_one();
    } else if (e == 1) {
      C.assign(A);
    } else {
      QQ<T> B;
      B.assign(A);
      ZZ exp;
      exp = e;
		
      if (e < 0) {
	B.invert(A);
	exp = -exp;
      }
      C.assign_one();
      while (exp > 0) {
	if (IsOdd(exp))
	  multiply(C, C, B);
	exp >>= 1;
	if (exp > 0)
	  square(B, B);
      }
    }
  }

  /* overloaded operators */
  template <class T> QQ<T> operator - (const QQ<T>& A) {
    QQ<T> C = A;
    C.negate();
    return C;
  }

  template <class T> QQ<T> operator + (const QQ<T>& A, const QQ<T>& B) {
    QQ<T> C;
    add(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator + (const QQ<T>& A, const T& B) {
    QQ<T> C;
    add(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator + (const T& A, const QQ<T>& B) {
    QQ<T> C;
    add(C, B, A);
    return C;
  }

  template <class T> QQ<T> operator - (const QQ<T>& A, const QQ<T>& B) {
    QQ<T> C;
    subtract(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator - (const QQ<T>& A, const T& B) {
    QQ<T> C;
    subtract(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator - (const T& A, const QQ<T>& B) {
    QQ<T> C;
    subtract(C, B, A);
    C.negate();
    return C;
  }

  template <class T> QQ<T> operator * (const QQ<T>& A, const QQ<T>& B) {
    QQ<T> C;
    multiply(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator * (const QQ<T>& A, const T& B) {
    QQ<T> C;
    multiply(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator * (const T& A, const QQ<T>& B) {
    QQ<T> C;
    multiply(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator / (const QQ<T>& A, const QQ<T>& B) {
    QQ<T> C;
    divide(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator / (const QQ<T>& A, const T& B) {
    QQ<T> C;
    divide(C, A, B);
    return C;
  }

  template <class T> QQ<T> operator / (const T& A, const QQ<T>& B) {
    QQ<T> C;
    divide(C, A, B);
    return C;
  }

  template <class T> const QQ<T>& operator += (QQ<T>& A, const QQ<T>& B) {
    add(A, A, B);
    return A;
  }

  template <class T> const QQ<T>& operator -= (QQ<T>& A, const QQ<T>& B) {
    subtract(A, A, B);
    return A;
  }

  template <class T> const QQ<T>& operator *= (QQ<T>& A, const QQ<T>& B) {
    multiply(A, A, B);
    return A;
  }

  template <class T> const QQ<T>& operator /= (QQ<T>& A, const QQ<T>& B) {
    divide(A, A, B);
    return A;
  }


  /* Input/Output */
  template <class T> std::istream& operator >> (std::istream& in, QQ<T>& X) {
    long n = 0;
    char c;

    in >> c;
    if (c != '(')
      ANTL_FATAL("ERROR:  QQ::operator>>::( expected");

    in >> c;
    while (c != ')' && n != 2) {
      in.putback(c);
      if (n == 0)
	in >> X.n;
      else
	in >> X.d;
      n++;
      in >> c;
      if (c == '/')
	in >> c;
    }
    if (::IsZero(X.d)) {
      ANTL_FATAL("ANTL: QQ divide by zero");
    }
	
    return in;
  }

  template <class T> std::ostream& operator << (std::ostream& out, const QQ<T>& X) {
    out << "(" << X.n << "/" << X.d << ")" << std::flush;
    return out;
  }



  template <class T> void QQ<T>::normalize() {
    // If the numerator is 0, always set denominator to 1
    if (NTL::IsZero(n)) {
      NTL::set(d);
      return;
    }
	
    // Make lowest terms
    T x;
    NTL::GCD(x, n, d%n);
    div(n,n,x);
    div(d,d,x);    
	
    // Make the denominator monic
    div(n,n,LeadCoeff(d));
    div(d,d,LeadCoeff(d));
  }

} // ANTL
