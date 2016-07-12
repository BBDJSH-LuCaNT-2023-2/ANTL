/**
 * @file xgcd_pseudo_impl.hpp
 * @author Michael Jacobson
 * @remark generic implementations of XGCD variations with pseudo-division
 */

//
// XGCD with only one inversion
//

template < class S, class T >
void XGCD_PSEUDO_work(T & G, T & X, T & Y, const T & A, const T & B)
{
  S d,z;

#ifdef DEBUG_PSEUDO
  cout << "\nIN XGCD_PSEUDO:  A = " << A << ", B = " << B << endl;
#endif

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
      PseudoDivRem<S,T>(d, q, u, u, v);
      swap(u, v);

      // update u sequence
      u0 = u2;

      // u2 = d u1 - q u2
      mul(temp, u2, q);
      mul(u2, u1, d);
      sub(u2, u2, temp);

      u1 = u0;


      // update v sequence
      v0 = v2;

      // v2 = d v1 - q v2
      mul(temp, v2, q);
      mul(v2, v1, d);
      sub(v2, v2, temp);

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

#ifdef DEBUG_PSEUDO
  if ((A*X + B*Y) != G) {
    cerr << "ERROR XGCD_PSEUDO" << endl;
    cerr << "AX + BY = " << A*X + B*Y << endl;
    cerr << "      G = " << G << endl;
    exit(1);
  }
#endif
}



//
// XGCD_LEFT versions with only one inversion
//

template < class S, class T >
void XGCD_LEFT_PSEUDO_work(T & G, T & X, const T & A, const T & B)
{
  S d,z;

#ifdef DEBUG_PSEUDO
  cout << "\nIN XGCD_LEFT_PSEUDO:  A = " << A << ", B = " << B << endl;
#endif

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
      PseudoDivRem<S,T>(d, q, u, u, v);
      swap(u, v);

      u0 = u2;

      // u2 = d u1 - q u2
      mul(temp, u2, q);
      mul(u2, u1, d);
      sub(u2, u2, temp);

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

#ifdef DEBUG_PSEUDO
  if ((A*X) % B != G % B) {
    cerr << "ERROR XGCD_LEFT_PSEUDO" << endl;
    cerr << "AX mod B = " << (A*X) % B << endl;
    cerr << "G % B = " << G % B << endl;
    exit(1);
  }
#endif
}


//
// Partial Euclidean algorithm (for NUCOMP) with pseudo-division
//  - Assumes that R1 is reduced modulo R2
//

template < class S, class T >
void XGCD_PARTIAL_PSEUDO_work(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag)
{
  S d,z;

  long e = deg(R2) + 1;
  T temp(INIT_SIZE, e), q(INIT_SIZE, e), r(INIT_SIZE,e),
    CC1(INIT_SIZE, e), CC2(INIT_SIZE, e);

  clear(CC2);
  set(CC1);
  NTL::negate(C1,C1);

  flag = false;
  set(z);

  while (deg (R1) > bound)
    {
      PseudoDivRem (d, q, R2, R2, R1);
      swap(R2,R1);
      mul(z,z,d);

      // r = d C2 - q C1
      mul(r,CC2,d);
      mul(temp,CC1,q);
      sub(r,r,temp);
      CC2 = CC1;
      CC1 = r;

      flag = !flag;
    }

  // remove inv(z)
  inv(z, z);
  mul(R1,R1,z);
  C2 = CC2;
  mul(C1,CC1,z);
}


//
// Partial Euclidean algorithm (for fast reduce) with pseudo-division
//  - Assumes that R1 is reduced modulo R2
//

template < class S, class T >
void XGCD_PARTIAL_REDUCE_PSEUDO_work(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even)
{
  S d,z;

  long e = deg(R2) + 1;
  T temp(INIT_SIZE, e), q(INIT_SIZE, e), r(INIT_SIZE,e),
    BB1(INIT_SIZE, e), BB2(INIT_SIZE, e);

  set(BB2);
  NTL::negate(BB2,BB2);
  clear(BB1);

  set(z);
  flag = false;

  while (deg (R1) > bound || (deg(R1) == bound && !even))
     {
      PseudoDivRem (d, q, R2, R2, R1);
      swap(R2,R1);
      mul(z,z,d);

      // r = d B2 - q B1
      mul(r,BB2,d);
      mul(temp,BB1,q);
      sub(r,r,temp);
      BB2 = BB1;
      BB1 = r;

      flag = !flag;
    }

  // remove inv(z)
  inv(z, z);
  mul(R1,R1,z);
  B2 = BB2;
  mul(B1,BB1,z);
}
