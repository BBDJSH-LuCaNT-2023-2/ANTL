/**
 * @file xgcd_plain_impl.hpp
 * @author Michael Jacobson
 * @remark generic implementations of XGCD_PARTIAL_PLAIN and XGCD_PARTIAL_REDUCE_PLAIN
 */

//
// XGCD_PLAIN
//

template < class S, class T >
void XGCD_PLAIN_work(T & G, T & X, T & Y, const T & A, const T & B)
{
#ifdef TRACE_XGCD
  cout << "--> IN XGCD_PLAIN" << endl;
#endif
  S z;

  if (IsZero(B)) {
    set(X);
    clear(Y);
    G = A;
  }
  else if (IsZero(A)) {
    clear(X);
    set(Y);
    G = B;
  }
  else {
    long e = max(deg(A), deg(B)) + 1;

    T temp(INIT_SIZE, e), u(INIT_SIZE, e), v(INIT_SIZE, e),
      u0(INIT_SIZE, e), u1(INIT_SIZE, e), u2(INIT_SIZE, e),
      v0(INIT_SIZE, e), v1(INIT_SIZE, e), v2(INIT_SIZE, e),
      q(INIT_SIZE, e);

    set(u1);
    clear(u2);
    clear(v1);
    set(v2);
    u = A; v = B;

    do {
      DivRem(q, u, u, v);
      swap(u, v);
      u0 = u2;
      mul(temp, q, u2);
      sub(u2, u1, temp);
      u1 = u0;
      v0 = v2;
      mul(temp, q, v2);
      sub(v2, v1, temp);
      v1 = v0;
    } while (!IsZero(v));

    G = u;
    X = u1;
    Y = v1;
  }

  if (IsZero(G)) return;
  if (IsOne(LeadCoeff(G))) return;

  /* make gcd monic */
  inv(z, LeadCoeff(G));
  mul(G, G, z);
  mul(X, X, z);
  mul(Y, Y, z);
}



//
// XGCD_LEFT_PLAIN
//

template < class S, class T >
void XGCD_LEFT_PLAIN_work(T & G, T & X, const T & A, const T & B)
{
#ifdef TRACE_XGCD
  cout << "--> IN XGCD_LEFT_PLAIN" << endl;
#endif
  S z;

  if (IsZero(B)) {
    set(X);
    G = A;
  }
  else if (IsZero(A)) {
    clear(X);
    G = B;
  }
  else {
    long e = max(deg(A), deg(B)) + 1;

    T temp(INIT_SIZE, e), u(INIT_SIZE, e), v(INIT_SIZE, e),
      u0(INIT_SIZE, e), u1(INIT_SIZE, e), u2(INIT_SIZE, e),
      q(INIT_SIZE, e);


    set(u1);
    clear(u2);
    u = A; v = B;

    do {
      DivRem(q, u, u, v);
      swap(u, v);
      u0 = u2;
      mul(temp, q, u2);
      sub(u2, u1, temp);
      u1 = u0;
    } while (!IsZero(v));

    G = u;
    X = u1;
  }

  if (IsZero(G)) return;
  if (IsOne(LeadCoeff(G))) return;

  /* make gcd monic */
  inv(z, LeadCoeff(G));
  mul(G, G, z);
  mul(X, X, z);
}




//
// Partial Euclidean algorithm (for NUCOMP)
//  - Assumes that R1 is reduced modulo R2
//

template < class T >
void XGCD_PARTIAL_PLAIN(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag)
{
  long e = deg(R2) + 1;
  T q(INIT_SIZE, e), r(INIT_SIZE,e);

  clear(C2);
  set(C1);
  NTL::negate(C1,C1);

  flag = false;

  while (deg (R1) > bound)
    {
      DivRem (q, R2, R2, R1);
      swap(R2,R1);

      // r = C2 - q C1
      mul(r,q,C1);
      sub(r,C2,r);
      C2 = C1;
      C1 = r;

      flag = !flag;
    }
}



//
// Partial Euclidean algorithm (for fast reduce)
//  - Assumes that R1 is reduced modulo R2
//

template < class T >
void XGCD_PARTIAL_REDUCE_PLAIN(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even)
{
  long e = deg(R2) + 1;
  T q(INIT_SIZE, e), r(INIT_SIZE,e);

  set(B2);
  NTL::negate(B2,B2);
  clear(B1);

  flag = false;

  while (deg (R1) > bound || (deg(R1) == bound && !even))
    {
      DivRem (q, R2, R2, R1);
      swap(R2,R1);

      // r = B2 - q B1
      mul(r,q,B1);
      sub(r,B2,r);
      B2 = B1;
      B1 = r;

      flag = !flag;
    }
}
