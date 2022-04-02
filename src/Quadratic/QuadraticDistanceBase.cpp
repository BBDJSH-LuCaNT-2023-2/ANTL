/**
 * @file QuadraticDistanceBase.cpp
 * @author Michael Jacobson
 * @version $Header$
 */

#include <ANTL/Quadratic/QuadraticDistanceBase.hpp>

using namespace std;

namespace ANTL {

// approximations of sqrt(Delta)
ZZ QuadraticDistanceBase::rootDs;
ZZ QuadraticDistanceBase::rootDp;
ZZ QuadraticDistanceBase::rootDpm;

ZZ QuadraticDistanceBase::xi;
ZZ QuadraticDistanceBase::yi;

// precision constants
long QuadraticDistanceBase::s;
long QuadraticDistanceBase::m;
long QuadraticDistanceBase::p;

// 2^p and log(2)
ZZ QuadraticDistanceBase::p2;
ZZ QuadraticDistanceBase::p21;
ZZ QuadraticDistanceBase::p22;
RR QuadraticDistanceBase::log2;

// output mode
const long QuadraticDistanceBase::REAL;
const long QuadraticDistanceBase::INTEGER;
long QuadraticDistanceBase::output_mode = 0;

//
// Constructor
//

QuadraticDistanceBase::QuadraticDistanceBase(const ZZ &nd, const ZZ &nk) {
  d = nd;
  k = nk;
}

//
// Copy Constructor
//

QuadraticDistanceBase::QuadraticDistanceBase(const QuadraticDistanceBase &qd) {
  d = qd.d;
  k = qd.k;
}

//
// set_precision
//
// Sets the precision variables using p such that 2^p > 41 S log_2 S, assuming
// that S = sqrt(Delta).
//

void QuadraticDistanceBase::set_precision(long Delta) {
  set_precision(to<ZZ>(Delta));
}

void QuadraticDistanceBase::set_precision(long long Delta) {
  set_precision(to<ZZ>(Delta));
}

void QuadraticDistanceBase::set_precision(const ZZ &Delta) {
  if (Delta < 0)
    return;

  ZZ S;

  SqrRoot(S, Delta);
  set_precision(Delta, S);
}

//
// set_precision
//
// Sets the precision variables using p such that 2^p > 41 S log_2 S
//

void QuadraticDistanceBase::set_precision(long Delta, const ZZ &S) {
  set_precision(to<ZZ>(Delta), S);
}

void QuadraticDistanceBase::set_precision(long long Delta, const ZZ &S) {
  set_precision(to<ZZ>(Delta), S);
}

void QuadraticDistanceBase::set_precision(const ZZ &Delta, const ZZ &S) {
  if (Delta < 0)
    return;

  // set RR precision
  if (NumBits(Delta) > 32) {
    RR::SetPrecision(NumBits(Delta));
    RR::SetOutputPrecision(NumBits(Delta) / 3);
  }
  log2 = log(to_RR(2));

  ZZ temp, C, rootD;
  long rootD_i;

  // 2^p > 46 B^2 log_2(B)
  temp = NumBits(S);
  if (temp > 16)
    temp = 16;
  p = NumBits(46 * temp * S * S);

  if (p < 32)
    p = 32;

  //
  // constants for floor (continued fraction expansion of sqrt(D))
  //

  // C = ceil(sqrt(2Delta (2^p Delta + 2 sqrt(Delta) + 2)))
  rootD = SqrRoot(Delta);
  temp = (Delta << p) + (rootD << 1) + 2;
  temp *= (Delta << 1);
  C = SqrRoot(temp) + 1;

  // compute convergents of sqrt(D)
  ZZ P, Q, NP, NQ, A0, A1, A2, B0, B1, B2, q;

  rootD_i = 0;
  clear(A2);
  set(A1);
  set(B2);
  clear(B1);

  set(Q);
  clear(P);
  do {
    if (Q > 0)
      q = (P + rootD) / Q;
    else
      q = (P + rootD + 1) / Q;

    NP = q * Q - P;
    NQ = (Delta - NP * NP) / Q;

    ++rootD_i;

    A0 = A2 + q * A1;
    B0 = B2 + q * B1;

    A2 = A1;
    A1 = A0;
    B2 = B1;
    B1 = B0;

    P = NP;
    Q = NQ;
  } while (B0 <= C);

  xi = A0;
  yi = B0;

  // p2 = 2^p
  set(p2);
  p2 = to<ZZ>(1) << p;

  // p21 = 2^(p+1)
  set(p21);
  p21 <<= (p + 1);

  // p22 = 2^(2 prec)
  set(p22);
  p22 <<= (p + p);

  // m = ceil(log2(2 D^1/2))
  m = NumBits(rootD << 1);

  // s = ceil(log2(D^sfrac))
  s = p + 1;

  // rootDs = floor(2^s sqrt(D))
  rootDs = (xi << s) / yi;

  // rootDp = floor(2^p sqrt(D))
  rootDp = (xi << p) / yi;

  // rootDpm = floor(2^(p+m) sqrt(D))
  rootDpm = (xi << (p + m)) / yi;
}

//
// set_output_mode
//
// Sets the mode used for output of QuadraticDistanceBase's
//    REAL (default) - output natural log of approx
//    INTEGER (needed for input) - output (d, k)
//

void QuadraticDistanceBase::set_output_mode(long mode) {
  if (mode == INTEGER)
    output_mode = INTEGER;
  else
    output_mode = REAL;
}

//
// assign_zero
//

void QuadraticDistanceBase::assign_zero() {
  clear(k);
  clear(d);
}

//
// assign_one
//

void QuadraticDistanceBase::assign_one() {
  // 1 = 2^0 2^p / 2^p --> k = 0, d = 2^p
  clear(k);
  d = p2;
}

//
// assign
//
/*
  void
  QuadraticDistanceBase::assign (const ZZ & d2, const ZZ & k2)
  {
  d = d2;
  k = k2;
  }
*/

//
// clear
//

void clear(QuadraticDistanceBase &qd) {
  clear(qd.k);
  clear(qd.d);
}

//
// set
//

void set(QuadraticDistanceBase &qd) {
  // 1 = 2^0 2^p / 2^p --> k = 0, d = 2^p
  clear(qd.k);
  qd.d = QuadraticDistanceBase::p2;
}

//
// assign(QuadraticDistanceBase)
//

void QuadraticDistanceBase::assign(const QuadraticDistanceBase &qd) {
  d = qd.d;
  k = qd.k;
}

//
// operator =
//

QuadraticDistanceBase &
QuadraticDistanceBase::operator=(const QuadraticDistanceBase &qd) {
  d = qd.d;
  k = qd.k;
  return *this;
}

//
// operator ==
//

bool operator==(const QuadraticDistanceBase &A,
                const QuadraticDistanceBase &B) {
  return A.d == B.d && A.k == B.k;
}

//
// operator !=
//

bool operator!=(const QuadraticDistanceBase &A,
                const QuadraticDistanceBase &B) {
  return A.d != B.d || A.k != B.k;
}

//
// swap
//

void swap(QuadraticDistanceBase &A, QuadraticDistanceBase &B) {
  QuadraticDistanceBase C;

  C = A;
  A = B;
  B = C;
}

//
// normalize
//
// Factors out powers of 2
//

void QuadraticDistanceBase::normalize() {
  ZZ temp = p2;
  long r = 0;

  while (temp < d) {
    ++r;
    temp <<= 1;
  }

  --r;

  temp = (to<ZZ>(1) << r) - 1;
  if ((d & temp) != 0) {
    d >>= r;
    ++d;
  } else
    d >>= r;

  k += r;
}

//
// invert(ZZ)
//
// Computes d^-1
//

void QuadraticDistanceBase::invert() {
  k = -k;
  d = (to<ZZ>(1) << (p << 1)) / d;
  ++d;
  normalize();
}

//
// multiply(ZZ)
//
// Computes d = ceil (d x)
//

void QuadraticDistanceBase::multiply(const ZZ &x) {
  d *= x;
  normalize();
}

//
// multiply(long)
//
// Computes d = ceil (d x)
//

void QuadraticDistanceBase::multiply(long x) {
  d *= to<ZZ>(x);
  normalize();
}

//
// multiply(long long)
//
// Computes d = ceil (d x)
//

void QuadraticDistanceBase::multiply(long long x) {
  d *= to<ZZ>(x);
  normalize();
}

//
// multiply_reduce(ZZ, ZZ)
//
// Computes d = ceil (d x / y)
//

void QuadraticDistanceBase::multiply_reduce(const ZZ &x, const ZZ &y) {
  d *= x;
  d /= y;
  ++d;
  normalize();
}

void QuadraticDistanceBase::multiply_reduce(long x, long y) {
  d *= to<ZZ>(x);
  d /= to<ZZ>(y);
  ++d;
  normalize();
}

void QuadraticDistanceBase::multiply_reduce(long long x, long long y) {
  d *= to<ZZ>(x);
  d /= to<ZZ>(y);
  ++d;
  normalize();
}

//
// multiply_rho(ZZ, ZZ)
//
// Computes d = ceil (d x / y)
//

void QuadraticDistanceBase::multiply_rho(const ZZ &x, const ZZ &y) {
  ZZ dmod = (rootDp + (x << p)) / y;

  ++dmod;
  d *= dmod;
  d >>= p;
  ++d;
  normalize();
}

void QuadraticDistanceBase::multiply_rho(long x, long y) {
  ZZ dmod = (rootDp + (to<ZZ>(x) << p)) / to<ZZ>(y);

  ++dmod;
  d *= dmod;
  d >>= p;
  ++d;
  normalize();
}

void QuadraticDistanceBase::multiply_rho(long long x, long long y) {
  ZZ dmod = (rootDp + (to<ZZ>(x) << p)) / to<ZZ>(y);

  ++dmod;
  d *= dmod;
  d >>= p;
  ++d;
  normalize();
}

//
// multiply_inverse_rho(ZZ, ZZ)
//
// Computes d = ceil (d x / y)
//

void QuadraticDistanceBase::multiply_inverse_rho(const ZZ &q, const ZZ &x,
                                                 const ZZ &y) {
  ZZ temp, dmod;
  long r = 0;

  set(temp);
  while (temp < q) {
    ++r;
    temp <<= 1;
  }

  dmod = (rootDpm - (x << (p + m))) / (y << (m - r));
  ++dmod;
  d *= dmod;

  set(temp);
  temp <<= ((p << 1) + 1);
  if (d > temp) {
    d >>= (p + 1);
    k -= (r - 1);
  } else {
    d >>= p;
    k -= r;
  }

  d = abs(d); // MDV if d is negative, the next loop never ends!
              //  Of course, I'm not sure this is 'correct'.
  while (d < p2) {
    d <<= 1;
    --k;
  }

  ++d;
}

void QuadraticDistanceBase::multiply_inverse_rho(long q, long x, long y) {
  ZZ temp, dmod;
  long r = 0;

  set(temp);
  while (temp < q) {
    ++r;
    temp <<= 1;
  }

  dmod = (rootDpm - (to<ZZ>(x) << (p + m))) / (to<ZZ>(y) << (m - r));
  ++dmod;
  d *= dmod;

  set(temp);
  temp <<= ((p << 1) + 1);
  if (d > temp) {
    d >>= (p + 1);
    k -= (r - 1);
  } else {
    d >>= p;
    k -= r;
  }

  d = abs(d);

  while (d < p2) {
    d <<= 1;
    --k;
  }

  ++d;
}

void QuadraticDistanceBase::multiply_inverse_rho(long long q, long long x,
                                                 long long y) {
  ZZ temp, dmod;
  long r = 0;

  set(temp);
  while (temp < to<ZZ>(q)) {
    ++r;
    temp <<= 1;
  }

  dmod = (rootDpm - (to<ZZ>(x) << (p + m))) / (to<ZZ>(y) << (m - r));
  ++dmod;
  d *= dmod;

  set(temp);
  temp <<= ((p << 1) + 1);
  if (d > temp) {
    d >>= (p + 1);
    k -= (r - 1);
  } else {
    d >>= p;
    k -= r;
  }

  d = abs(d);

  while (d < p2) {
    d <<= 1;
    --k;
  }

  ++d;
}

//
// divide(ZZ)
//
// Computes d = ceil (d / y)
//

void QuadraticDistanceBase::divide(const ZZ &x) {
  ZZ temp;

  temp = d / x;
  while (temp < p2) {
    d <<= 1;
    --k;
    temp = d / x;
  }

  d = temp + 1;
  normalize();
}

void QuadraticDistanceBase::divide(long x) {
  ZZ temp;

  temp = d / x;
  while (temp < p2) {
    d <<= 1;
    --k;
    temp = d / to<ZZ>(x);
  }

  d = temp + 1;
  normalize();
}

void QuadraticDistanceBase::divide(long long x) {
  ZZ temp, xZZ;

  xZZ = to<ZZ>(x);
  temp = d / xZZ;
  while (temp < p2) {
    d <<= 1;
    --k;
    temp = d / xZZ;
  }

  d = temp + 1;
  normalize();
}

//
// divide_reduce()
//
// Task:
//      divides the nucomp distance approximation by the reduction distance
//

void QuadraticDistanceBase::reduce_div(const ZZ &x) {
  register long t;
  ZZ r, temp2;

  t = 0;
  temp2 = (d << 3);
  while (temp2 <= x) {
    temp2 <<= 1;
    ++t;
  }
  temp2 <<= p;
  DivRem(d, r, temp2, x);
  if (!IsZero(r))
    ++d;
  k -= t;

  normalize();
}

//
// multiply(QuadraticDistanceBase, QuadraticDistanceBase, ZZ, ZZ)
//
// Computes the product of the QuadraticDistanceBase with the
// QuadraticDistanceBase x/y
//

void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
              const ZZ &x, const ZZ &y) {
  ZZ dmod = (x << QuadraticDistanceBase::p) / y;

  ++dmod;
  C.d = abs(A.d * dmod);
  C.d >>= QuadraticDistanceBase::p;
  ++C.d;
  C.k = A.k;

#if defined(DEBUG_COMP) || defined(DEBUG_NUCOMP)
  cout << "\ndmod = " << dmod << ", dist = " << C.d << endl;
#endif
}

void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
              const long &x, const long &y) {
  ZZ dmod = (to<ZZ>(x) << QuadraticDistanceBase::p) / to<ZZ>(y);

  ++dmod;
  C.d = abs(A.d * dmod);
  C.d >>= QuadraticDistanceBase::p;
  ++C.d;
  C.k = A.k;

#if defined(DEBUG_COMP) || defined(DEBUG_NUCOMP)
  cout << "\ndmod = " << dmod << ", dist = " << C.d << endl;
#endif
}

void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
              const long long &x, const long long &y) {
  ZZ dmod = (to<ZZ>(x) << QuadraticDistanceBase::p) / to<ZZ>(y);

  ++dmod;
  C.d = abs(A.d * dmod);
  C.d >>= QuadraticDistanceBase::p;
  ++C.d;
  C.k = A.k;

#if defined(DEBUG_COMP) || defined(DEBUG_NUCOMP)
  cout << "\ndmod = " << dmod << ", dist = " << C.d << endl;
#endif
}

//
// multiply(QuadraticDistanceBase, QuadraticDistanceBase, qo_distance)
//
// Computes the product of two QuadraticDistanceBases
//

void multiply(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
              const QuadraticDistanceBase &B) {

  /*
    C.d = (A.d * B.d) >> QuadraticDistanceBase::p;
    C.k = A.k + B.k;
    C.normalize();
  */

  ZZ dmod, e, ee;

  dmod = (A.d * B.d) << 1;
  e = (dmod + QuadraticDistanceBase::p2) >> (QuadraticDistanceBase::p + 1);
  ee = (dmod + QuadraticDistanceBase::p21) >> (QuadraticDistanceBase::p + 2);
  C.k = A.k + B.k;
  if (e <= QuadraticDistanceBase::p21)
    C.d = e;
  else if (ee >= (QuadraticDistanceBase::p2 + 1)) {
    C.d = ee;
    ++C.k;
  } else
    C.d = QuadraticDistanceBase::p21;
}

//
// square(QuadraticDistanceBase, qo_distace)
//
// Computes the square of the QuadraticDistanceBase
//

void square(QuadraticDistanceBase &C, const QuadraticDistanceBase &A) {

  /*
    C.d = (A.d * A.d) >> QuadraticDistanceBase::p;
    C.k = (A.k << 1);
    C.normalize();
  */

  ZZ dmod, e, ee;

  dmod = (A.d * A.d) << 1;
  e = (dmod + QuadraticDistanceBase::p2) >> (QuadraticDistanceBase::p + 1);
  ee = (dmod + QuadraticDistanceBase::p21) >> (QuadraticDistanceBase::p + 2);
  C.k = (A.k << 1);
  if (e <= QuadraticDistanceBase::p21)
    C.d = e;
  else if (ee >= (QuadraticDistanceBase::p2 + 1)) {
    C.d = ee;
    ++C.k;
  } else
    C.d = QuadraticDistanceBase::p21;
}

//
// divide(QuadraticDistanceBase, QuadraticDistanceBase, QuadraticDistanceBase)
//
// Computes the product of two QuadraticDistanceBases
//

void divide(QuadraticDistanceBase &C, const QuadraticDistanceBase &A,
            const QuadraticDistanceBase &B) {
  C.d = (A.d << QuadraticDistanceBase::p) / B.d;
  C.k = A.k - B.k;
  C.normalize();
}

//
// operator >>
//

std::istream &operator>>(std::istream &in, QuadraticDistanceBase &A) {
  long n = 0;
  char c;
  ZZ ibuf[2];

  in >> c;
  if (c != '(') {
    cerr << "ERROR:  QuadraticDistanceBase::operator>>::( expected" << endl;
    exit(1);
  }

  in >> c;
  while (c != ')' && n != 2) {
    in.putback(c);
    in >> ibuf[n];
    n++;
    in >> c;
    if (c == ',')
      in >> c;
  }

  A.d = ibuf[0];
  A.k = ibuf[1];

  return in;
}

//
// operator <<
//

std::ostream &operator<<(std::ostream &out, const QuadraticDistanceBase &A) {
  if (QuadraticDistanceBase::output_mode == QuadraticDistanceBase::INTEGER)
    out << "(" << A.d << ", " << A.k << ")" << flush;
  else
    out << "(" << A.d << ", " << A.k << ", " << A.get_log() << ")" << flush;

  return out;
}

} // namespace ANTL
