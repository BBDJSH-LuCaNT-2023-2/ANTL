/**
 * @file qo_reduce_fast_ZZ.cpp
 * @author Michael Jacobson
 * @remark Specialization of the qo_reduce_fast class for ZZ. 
 */

#include <Quadratic/Reduce/qo_reduce_fast.hpp>
#include <NTL/RR.h>


template <>
void 
qo_reduce_fast<ZZ>::init(const ZZ & Din, const ZZ & hin, long gin)
{
  qo_reduce<ZZ>::init(Din,hin,0);
  SQRT_DELTA = FloorToZZ(sqrt(abs(to_RR(Delta))));
}


template <>
void
qo_reduce_fast<ZZ>::reduce(QuadraticIdealBase<ZZ> & A)
{
  static ZZ a, b, c, na, nb, q, r, a2, temp;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();

  if (a < SQRT_DELTA) {
    // just normalize - close to reduced
    if (b <= -a || b > a) {
      LeftShift(a2,a,1);
  
      // q = b/2a
      DivRem(q, r, b, a2);

      if (r > a) {
	sub(r,r,a2);
	++q;
      }

      // c -= q (b + r) / 2
      add(temp,b,r);
      RightShift(temp,temp,1);
      mul(temp,temp,q);
      sub(c,c,temp);

      // b = r
      b = r;
    }
  }
  else {
    // reduce
    ZZ RR,R,CC,C,N,oa;

    LeftShift(oa,abs(a),1);

    RR = oa;
    R = b;

    // use bound sqrt(a) D^1/4 / 2
    mul(N, SQRT_DELTA, a);
    RightShift(N,N,1);
    SqrRoot(N, N);

    XGCD_PARTIAL(RR,R,CC,C,N);

    // a = (-1)^(i+1) (R^2 - Delta C^2) / (4 a0)
    sqr(a,R);
    sqr(temp,C);
    mul(temp,temp,Delta);
    sub(a,a,temp);
    div(a,a,oa);
    RightShift(a,a,1);
    if (C < 0)
      NTL::negate(a,a);

    // b = (R + a BB) / B
    mul(temp,a,CC);
    LeftShift(temp,temp,1);
    add(b,R,temp);
    div(b,b,C);

    if (a < 0)
      NTL::negate(a,a);
    LeftShift(oa,a,1);
    rem(b,b,oa);
    if (b > a)
      sub(b,b,oa);

    sqr(c,b);
    sub(c,c,Delta);
    div(c,c,a);
    RightShift(c,c,2); 
  }

  // one additional reduction step if necessary
  while (a > c) {
    na = c;

    LeftShift(a2,na,1);

    // -b = 2q * na + nb
    NTL::negate(temp,b);
    DivRem (q, nb, temp, a2);

    if (nb > na)
      {
	sub(nb,nb,a2);
	++q;
      }

    // c = a - q * (nb - b)/2
    sub(temp,nb,b);
    RightShift(temp,temp,1);
    mul(temp,temp,q);
    sub(c,a,temp);

    b = nb;
    a = na;
  }

  // account for special case
  if ((a == c) && (b < 0))
    NTL::negate(b,b);

  A.assign(a,b,c);
}
