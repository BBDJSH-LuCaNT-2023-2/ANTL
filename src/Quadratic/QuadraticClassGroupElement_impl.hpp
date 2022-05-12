/**
 * @file QuadraticClassGroupElement_impl.hpp
 * @author Michael Jacobson
 * @remarks This file is to be included from QuadraticClassGroupElement.h only.
 * @version $Header$
 */

namespace ANTL {

//
// Public functions
//

//
// constructor
//

template <class T>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement()
    : QuadraticIdealBase<T>() {}

//
// constructor:
//    initialize with a reference to a Quadratic Order.
//

template <class T>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement(
    QuadraticOrder<T> &qo_a) {
  QO = &qo_a;
}


//
// constructor:
//    initialize with a copy of A.
//

template <class T>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement(
    const QuadraticIdealBase<T> &A)
    : QuadraticIdealBase<T>(A) {
  this->reduce();
}

//
// constructor:
//    initialize with a HashEntry.
//

template <class T>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement(const HashEntry<T> &A)
    : QuadraticIdealBase<T>(A) {
  this->reduce();
}

//
// constructor:
//    initialize with a HashEntry.
//

template <class T>
template <class S>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement(
    const HashEntryInt<T, S> &A)
    : QuadraticIdealBase<T>(A) {
  this->reduce();
}

//
// constructor:
//    initialize with a HashEntry.
//

template <class T>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement(
    const HashEntryVec<T> &A)
    : QuadraticIdealBase<T>(A) {
  this->reduce();
}

//
// constructor:
//    initialize with a copy of A.
//

template <class T>
QuadraticClassGroupElement<T>::QuadraticClassGroupElement(
    const QuadraticClassGroupElement<T> &A) {
  this->a = A.a;
  this->b = A.b;
  this->c = A.c;
}

//
// destructor
//

template <class T>
QuadraticClassGroupElement<T>::~QuadraticClassGroupElement() {}

/*
//
// QuadraticClassGroupElement<T>::set_current_order()
//
// Task:
//      sets the current quadratic order and the other static variables
//

template < class T >
void
QuadraticClassGroupElement<T>::set_current_order(const T & newDelta)
{
QuadraticIdealBase<T>::set_current_order(newDelta);
}


template < class T >
void
QuadraticClassGroupElement<T>::set_current_order(const T & newDelta, const T &
newh)
{
QuadraticIdealBase<T>::set_current_order(newDelta,newh);
}
*/

//
// assign_prime()
//
// Task:
//      computes a reduced representative of the equivalence class containing
//      the ideal lying over the prime p.  If such an ideal doesnot exist,
//      false is returned.
//

template <class T>
bool QuadraticClassGroupElement<T>::assign_prime(const T &p) {
  bool OK = QuadraticIdealBase<T>::assign_prime(p);

  if (OK)
    this->reduce();

  return OK;
}

//
// QuadraticClassGroupElement<T>::assign(T, T)
//
// Task:
//      set to the form with the given coeffs.  Returns false if not possible.
//

template <class T>
bool QuadraticClassGroupElement<T>::assign(const T &a2, const T &b2) {
  bool OK = QuadraticIdealBase<T>::assign(a2, b2);

  if (OK)
    this->reduce();

  return OK;
}

//
// QuadraticClassGroupElement<T>::assign(HashEntry<T>)
//
// Task:
//      set to a copy of B.
//

template <class T>
void QuadraticClassGroupElement<T>::assign(const HashEntry<T> &B) {
  QuadraticIdealBase<T>::assign(B);
  this->reduce();
}

//
// QuadraticClassGroupElement<T>::assign(HashEntryInt<T>)
//
// Task:
//      set to a copy of B.
//

template <class T>
template <class S>
void QuadraticClassGroupElement<T>::assign(const HashEntryInt<T, S> &B) {
  std::cout << "assign: assigning HashEntryInt<T, S>!" << std::endl;
  QuadraticIdealBase<T>::assign(B);
  std::cout << "assign: reducing!" << std::endl;
  this->reduce();
}

//
// QuadraticClassGroupElement<T>::assign(HashEntryVec<T>)
//
// Task:
//      set to a copy of B.
//

template <class T>
void QuadraticClassGroupElement<T>::assign(const HashEntryVec<T> &B) {
  QuadraticIdealBase<T>::assign(B);
  this->reduce();
}

//
// QuadraticClassGroupElement<T>::assign(QuadraticClassGroupElement<T>)
//
// Task:
//      set to a copy of B.
//

template <class T>
void QuadraticClassGroupElement<T>::assign(
    const QuadraticClassGroupElement<T> &B) {
  this->a = B.a;
  this->b = B.b;
  this->c = B.c;
}


//
// operator =
//
// Task:
//      make a copy of an existing QuadraticClassGroupElement<T>
//

template <class T>
QuadraticClassGroupElement<T> &QuadraticClassGroupElement<T>::operator=(
    const QuadraticClassGroupElement<T> &A) {
  assign(A);
  return *this;
}

template <class T>
QuadraticClassGroupElement<T> conjugate(const QuadraticClassGroupElement<T> &C){
  //std::cout << "Conjugating!" << std::endl;
  QuadraticClassGroupElement<T> conj_C(C);

  conj_C.a = C.a;
  conj_C.b = -C.b;
  conj_C.c = C.c;

  conj_C.reduce();
  return conj_C;
}

//
// multiply_imag()
//
// Task:
//      multiplies two real ideal equivalence classes
//

template <class T>
void multiply_imag(QuadraticClassGroupElement<T> &C,
                   const QuadraticClassGroupElement<T> &A,
                   const QuadraticClassGroupElement<T> &B) {
  T junk;
  junk = multiply_base(C, A, B);
  C.reduce_imag_mult();

#ifdef DEBUG_COMP
  cout << "COMPOSITE:  " << C << endl;
  C.test_form((char *)"MULTIPLY_IMAG");
#endif
}

//
// multiply_real()
//
// Task:
//      multiplies two real ideal equivalence classes
//

template <class T>
void multiply_real(QuadraticClassGroupElement<T> &C,
                   const QuadraticClassGroupElement<T> &A,
                   const QuadraticClassGroupElement<T> &B) {
  T junk = multiply_base(C, A, B);
  C.reduce_real();

#ifdef DEBUG_COMP
  cout << "COMPOSITE:  " << C << endl;
  C.test_form((char *)"MULTIPLY_REAL");
#endif
}

//
// multiply()
//
// Task:
//      multiplies two real ideal equivalence classes
//

template <class T>
void mul(QuadraticClassGroupElement<T> &C,
              const QuadraticClassGroupElement<T> &A,
              const QuadraticClassGroupElement<T> &B) {
  C.QO->get_mul_best()->multiply(C, A, B);
  C.reduce();

#ifdef DEBUG_COMP
  cout << "COMPOSITE:  " << C << endl;
  C.test_form((char *)"MULTIPLY");
#endif
}

//
// square_imag()
//
// Task:
//      squares an ideal equivalence class
//

template <class T>
void square_imag(QuadraticClassGroupElement<T> &C,
                 const QuadraticClassGroupElement<T> &A) {
  T junk;
  junk = square_base(C, A);
  C.reduce_imag_mult();

#ifdef DEBUG_COMP
  cout << "COMPOSITE:  " << C << endl;
  C.test_form((char *)"SQUARE_IMAG");
#endif
}

//
// square_real()
//
// Task:
//      squares an ideal equivalence class
//

template <class T>
void square_real(QuadraticClassGroupElement<T> &C,
                 const QuadraticClassGroupElement<T> &A) {
  T junk = square_base(C, A);
  C.reduce_real();

#ifdef DEBUG_COMP
  cout << "COMPOSITE:  " << C << endl;
  C.test_form((char *)"SQUARE_REAL");
#endif
}

//
// square()
//
// Task:
//      squares an ideal equivalence class
//

template <class T>
void square(QuadraticClassGroupElement<T> &C,
            const QuadraticClassGroupElement<T> &A) {
  C.QO->get_sqr_best()->square(C, A);
  C.reduce();

#ifdef DEBUG_COMP
  cout << "COMPOSITE:  " << C << endl;
  C.test_form((char *)"SQUARE");
#endif
}

//
// power_imag(T)
//
// Task:
//      computes A^i using binary exponentiation
//

template <class T>
void power_imag(QuadraticClassGroupElement<T> &C,
                const QuadraticClassGroupElement<T> &A, const ZZ &n) {
  QuadraticClassGroupElement<T> B;
  long i, k = 0;
  ZZ j, ex;

  if (IsZero(n)) {
    C.assign_one();
    return;
  }

  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  if (n > 0) {
    B.assign(A);
    C.assign(A);
  } else {
    B = conjugate(A);
    C = conjugate(A);
  }

  for (i = 1; i <= k; ++i) {
    square_imag(C, C);
    if (IsOdd(j))
      multiply_imag(C, C, B);
    j >>= 1;
  }
}

//
// power_real(T)
//
// Task:
//      computes A^i using binary exponentiation
//

template <class T>
void power_real(QuadraticClassGroupElement<T> &C,
                const QuadraticClassGroupElement<T> &A, const ZZ &n) {
  QuadraticClassGroupElement<T> B;
  long i, k = 0;
  ZZ j, ex;

  if (IsZero(n)) {
    C.assign_one();
    return;
  }

  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  if (n > 0) {
    B.assign(A);
    C.assign(A);
  } else {
    B = conjugate(A);
    C = conjugate(A);
  }

  for (i = 1; i <= k; ++i) {
    square_real(C, C);
    if (IsOdd(j))
      multiply_real(C, C, B);
    j >>= 1;
  }
}

//
// power(T)
//
// Task:
//      computes A^i using binary exponentiation
//

template <class T>
void power(QuadraticClassGroupElement<T> &C,
           const QuadraticClassGroupElement<T> &A, const ZZ &n) {
  QuadraticClassGroupElement<T> B;
  long i, k = 0;
  ZZ j, ex;

  if (IsZero(n)) {
    C.assign_one();
    return;
  }

  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  if (n > 0) {
    B.assign(A);
    C.assign(A);
  } else {
    B = conjugate(A);
    C = conjugate(A);
  }

  for (i = 1; i <= k; ++i) {
    square(C, C);
    if (IsOdd(j))
      multiply(C, C, B);
    j >>= 1;
  }
}

//
// nupower_imag(T)
//
// Task:
//      computes A^i using binary exponentiation
//

template <class T>
void nupower_imag(QuadraticClassGroupElement<T> &C,
                  const QuadraticClassGroupElement<T> &A, const ZZ &n) {
  QuadraticClassGroupElement<T> B;
  ZZ j, ex;
  long i, k = 0;

  if (IsZero(n)) {
    C.assign_one();
    return;
  }

  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  if (n > 0) {
    B.assign(A);
    C.assign(A);
  } else {
    B = conjugate(A);
    C = conjugate(A);
  }

  for (i = 1; i <= k; ++i) {
    nudupl_imag(C, C);
    if (IsOdd(j))
      nucomp_imag(C, C, B);
    j >>= 1;
  }
}

//
// nupower_real(T)
//
// Task:
//      computes A^i using binary exponentiation
//

template <class T>
void nupower_real(QuadraticClassGroupElement<T> &C,
                  const QuadraticClassGroupElement<T> &A, const ZZ &n) {
  QuadraticClassGroupElement<T> B;
  ZZ j, ex;
  long i, k = 0;

  if (IsZero(n)) {
    C.assign_one();
    return;
  }

  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  if (n > 0) {
    B.assign(A);
    C.assign(A);
  } else {
    B = conjugate(A);
    C = conjugate(A);
  }

  for (i = 1; i <= k; ++i) {
//     nudupl_real(C, C);
    mul(C, C, C);
    if (IsOdd(j))
//       nucomp_real(C, C, B);
      mul(C, C, B);
    j >>= 1;
  }
}

//
// nupower(T)
//
// Task:
//      computes A^i using binary exponentiation
//

template <class T>
void nupower(QuadraticClassGroupElement<T> &C,
             const QuadraticClassGroupElement<T> &A, const ZZ &n) {
  QuadraticClassGroupElement<T> B;
  ZZ j, ex;
  long i, k = 0;

  if (IsZero(n)) {
    C.assign_one();
    return;
  }

  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  if (n > 0) {
    B.assign(A);
    C.assign(A);
  } else {
    B = conjugate(A);
    C = conjugate(A);
  }

  for (i = 1; i <= k; ++i) {
    nudupl(C, C);
    if (IsOdd(j))
      nucomp(C, C, B);
    j >>= 1;
  }
}

//
//
// operator *
//
// Task:
//      multiplies A and B
//

template <class T>
QuadraticClassGroupElement<T>
operator*(const QuadraticClassGroupElement<T> &A,
          const QuadraticClassGroupElement<T> &B) {
  QuadraticClassGroupElement<T> C;

  nucomp(C, A, B);
  return C;
}

//
// operator *=
//
// Task:
//      *this = *this * A
//

template <class T>
QuadraticClassGroupElement<T> &QuadraticClassGroupElement<T>::operator*=(
    const QuadraticClassGroupElement<T> &A) {
  nucomp(*this, *this, A);
  return *this;
}

//
// swap
//
// Task:
//      swaps A and B
//

template <class T>
void swap(QuadraticClassGroupElement<T> &A, QuadraticClassGroupElement<T> &B) {
  QuadraticClassGroupElement<T> C;

  C.assign(A);
  A.assign(B);
  B.assign(C);
}

//
// operator >>
//
// Task:
//      inputs a QuadraticClassGroupElement<T> from the istream in.
//

template <class T>
std::istream &operator>>(std::istream &in, QuadraticClassGroupElement<T> &A) {
  long n = 0;
  char c;
  T ibuf[3];

  in >> c;
  if (c != '(') {
    cout << "ERROR:  QuadraticClassGroupElement<T>::operator>>::( expected"
         << endl;
    exit(1);
  }

  in >> c;
  while (c != ')' && n != 3) {
    in.putback(c);
    in >> ibuf[n];
    n++;
    in >> c;
    if (c == ',')
      in >> c;
  }

  A.a = ibuf[0];
  A.b = ibuf[1];
  A.c = ibuf[2];
  A.reduce();

  return in;
}

} // namespace ANTL
