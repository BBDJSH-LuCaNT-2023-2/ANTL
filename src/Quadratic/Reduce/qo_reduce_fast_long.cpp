/**
 * @file qo_reduce_fast_long.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_reduce_fast class for long. 
 */

#include <Quadratic/Reduce/qo_reduce_fast.hpp>
#include <NTL/RR.h>


template <>
void 
qo_reduce_fast<long>::init(const long & Din, const long & hin, long gin)
{
  qo_reduce<long>::init(Din,hin,0);
  SQRT_DELTA = FloorTolong(sqrt(abs(to_RR(Delta))));
}


template <>
void
qo_reduce_fast<long>::reduce(QuadraticIdealBase<long> & A)
{
  static long a, b, c, na, nb, q, r, a2, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  if (a < SQRT_DELTA) {
    // just normalize - close to reduced
    if (b <= -a || b > a) {
      a2 = a << 1;
  
      // q = b/2a
      DivRem(q, r, b, a2);

      if (r > a) {
	r -= a2;
	++q;
      }

      // c -= q (b + r) / 2
      c -= q*(b+r) >> 1;

      // b = r
      b = r;
    }
  }
  else {
    // reduce
    long RR,R,CC,C,N,oa;

    oa = abs(a) << 1;

    RR = oa;
    R = b;

    // use bound sqrt(a) D^1/4 / 2
    mul(N, SQRT_DELTA, a);
    RightShift(N,N,1);
    SqrRoot(N, N);

    XGCD_PARTIAL(RR,R,CC,C,N);

    // a = (-1)^(i+1) (R^2 - Delta C^2) / (4 a0)
    a = ((R*R - Delta*C*C) / a0) >> 1;   // MJJ:  ERROR?
    if (C < 0)
      a = -a;

    // b = (R + a BB) / B
    b = (R + (a*CC << 1))) / C;

    if (a < 0)
      a = -a;
    oa = a << 1;
    b %= oa;
    if (b > a)
      b -= oa;

    c = ((b*b - Delta) / a) >> 2;
  }

  // one additional reduction step if necessary
  while (a > c) {
    na = c;

    a2 = na << 1;

    // -b = 2q * na + nb
    DivRem (q, nb, -b, a2);

    if (nb > na)
      {
	nb -= a2;
	++q;
      }

    // c = a - q * (nb - b)/2
    c = a - q * (nb -b) >> 1;

    b = nb;
    a = na;
  }

  // account for special case
  if ((a == c) && (b < 0))
    b = -b;

  A.assign(a,b,c);
}
