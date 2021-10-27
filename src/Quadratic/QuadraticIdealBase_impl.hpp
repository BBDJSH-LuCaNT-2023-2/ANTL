/**
 * @file QuadraticIdealBase_impl.hpp
 * @author Michael Jacobson
 * @remarks Class representing primitive ideal.
 */

// constructor
template <class T> QuadraticIdealBase<T>::QuadraticIdealBase (QuadraticOrder<T> & inQO) {
  NTL::clear(a);
  NTL::clear(b);
  NTL::clear(c);
  QO = &inQO;
}

// destructor
template <class T> QuadraticIdealBase<T>::~QuadraticIdealBase () {}

// Getters and Setters
template <class T> T QuadraticIdealBase<T>::get_a () const {return a;}
template <class T> T QuadraticIdealBase<T>::get_b () const {return b;}
template <class T> T QuadraticIdealBase<T>::get_c () const {return c;}

template <class T> void QuadraticIdealBase<T>::set_a  (T x) {a = x;}
template <class T> void QuadraticIdealBase<T>::set_b  (T x) {b = x;}
template <class T> void QuadraticIdealBase<T>::set_c  (T x) {c = x;}

template <class T> QuadraticOrder<T>* QuadraticIdealBase<T>::get_QO () const {return QO;}
template <class T> void QuadraticIdealBase<T>::set_QO (QuadraticOrder<T> *x) {QO = x;}

// template <class T> void QuadraticIdealBase<T>::test_ideal (string msg) {
//   T temp = b * b - 4 * a * c;
//   if (temp != QO->getDiscriminant())
//     {
//       cout << "ERROR " << msg << "!  wrong getDiscriminant!" << endl;
//       cout << "a = " << a << ", b = " << b << ", c = " << c << endl;
//       cout << "Delta = " << QO->getDiscriminant() << endl;
//       cout << "b^2 - ac = " << temp << endl;
//       cout << "(b^2 - Delta) % a = " << (b*b - QO->getDiscriminant()) % a << endl;
//       exit (1);
//     }
// }

// QuadraticIdealBase<T>::assign_one()
//
// Task: set to the unit ideal of the quadratic_order
template <class T> void QuadraticIdealBase<T>::assign_one () {
  set (a);
  clear (b);
  c = -QO->getDiscriminant();
}

// assign_prime
//
// Task: computes a reduced representative of the equivalence class containing
//      the ideal lying over the prime p.  If such an ideal doesnot exist,
//      false is returned.
template <class T> bool QuadraticIdealBase<T>::assign_prime (const T & p) {
  T temp, Dp;
  long jac;

  if (!DetIrredTest (p))
    return false;

  jac = ressol (temp, QO->getDiscriminant() % p, p);

  if (jac < 0)
    return false;

  if (jac == 0)
    {
      temp = QO->getDiscriminant() % (p * p);
      if (IsZero (temp))
    return false;
      else
    clear (temp);
    }

  a = p;
  rem(b,temp,a);

  // c = (b^2 - Delta) / a
  sqr(c,b);
  sub(c,c,QO->getDiscriminant());
  div(c,c,a);

  return true;
}

// QuadraticIdealBase<T>::assign(T,T,T)
//
// Task: set to a copy of B.
template <class T> void QuadraticIdealBase<T>::assign (const T & na, const T & nb, const T & nc) {
  a = na;
  b = nb;
  c = nc;
}

// QuadraticIdealBase<T>::assign(QuadraticIdealBase<T>)
//
// Task: set to a copy of B.
template <class T> void QuadraticIdealBase<T>::assign (const QuadraticIdealBase<T> &B) {
  a = B.a;
  b = B.b;
  c = B.c;
}

// operator =
//
// Task: make a copy of an existing QuadraticIdealBase<T>
template <class T> QuadraticIdealBase<T> & QuadraticIdealBase<T>::operator = (const QuadraticIdealBase<T> &A) {
  assign (A);
  return *this;
}

// conjugate()
//
// Task:computes the conjugate of A
template <class T> void conjugate (QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A) {
  C.a = A.a;
  C.b = -A.b;
  C.c = A.c;
}

template <class T> void ANTL::mul(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B) {
  C.QO->get_mul_best()->multiply(C, A, B);
}

template <class T> void mul(QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B) {
  C.QO->get_mul_best()->multiply(C,gamma,A,B);
}

template <class T> void sqr(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A) {
  C.QO->get_sqr_best()->sqr(C,A);
}


template <class T> void sqr(QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A) {
  C.QO->get_sqr_best()->sqr(C,gamma,A);
}

template <class T> void cube(QuadraticIdealBase<T> &C, const QuadraticIdealBase<T> &A) {
  C.QO->get_cube_best()->cube(C,A);
}

template <class T>void cube(QuadraticIdealBase<T> &C, ANTL::QuadraticNumber<T> & gamma, const QuadraticIdealBase<T> &A) {
  C.QO->get_cube_best()->cube(C,gamma,A);
}

template <class T> void QuadraticIdealBase<T>::reduce() {
  QO->get_red_best()->reduce(*this);
}

template <class T> void QuadraticIdealBase<T>::reduce(ANTL::QuadraticNumber<T> & gamma) {
  QO->get_red_best()->reduce(*this,gamma);
}

// QuadraticIdealBase<T>::IsOne()
//
// Task: tests if the ideal is the unit ideal
template <class T> bool QuadraticIdealBase<T>::IsOne () const {
  return (IsOne (a));
}

// QuadraticIdealBase<T>::IsEqual()
//
// Task: tests if the ideals are equal
template <class T> bool QuadraticIdealBase<T>::IsEqual (const QuadraticIdealBase<T> &B) const {
  return ((a == B.a) && (b == B.b) && QO->IsEqual(*B.QO));
}

// operator ==
//
// Task: tests if A and B are equal
template < class T > bool operator == (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B) {
  return (A.IsEqual (B));
}

// operator !=
//
// Task: tests if A and B are not equal
template <class T> bool operator != (const QuadraticIdealBase<T> &A, const QuadraticIdealBase<T> &B) {
  return (!A.IsEqual (B));
}

// operator >>
//
// Task: inputs a QuadraticIdealBase<T> from the istream in.
template <class T> std::istream & operator >> (std::istream & in, QuadraticIdealBase<T> &A) {
  long n = 0;
  char c;
  T ibuf[3];

  in >> c;
  if (c != '(') {
      cout << "ERROR:  QuadraticIdealBase<T>::operator>>::( expected" << endl;
      exit (1);
    }

  in >> c;
  while (c != ')' && n != 3) {
      in.putback (c);
      in >> ibuf[n];
      n++;
      in >> c;
      if (c == ',')
    in >> c;
    }

  A.a = ibuf[0];
  A.b = ibuf[1];
  A.c = ibuf[2];

  return in;
}

// operator <<
//
// Task: outputs a QuadraticIdealBase<T> to the ostream out.
template <class T> std::ostream & operator << (std::ostream & out, const QuadraticIdealBase<T> &A) {
  out << "(" << A.a << ", " << A.b << ", " << A.c << ")";
  return out;
}
