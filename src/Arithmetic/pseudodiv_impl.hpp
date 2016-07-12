/**
 * @file pseudodiv_impl.hpp
 * @author Michael Jacobson
 * @remark generic implementations of pseudodivision algorithm
 */

//
// Pseudo division with remainder (for 1 inversion XGCD routines)
//

template < class S, class T >
void PseudoDivRem(S & d, T & q, T & r, const T & a, const T & b)
{
  // Calls PseudoDivRem_reduce by default.  Calls _reduce or _accumulate
  // version for ZZ_pX and GF2EX based on pre-computed threshold.

  PseudoDivRem_reduce(d,q,r,a,b);
}



template < class S, class T >
void PseudoDivRem_reduce(S & d, T & q, T & r, const T & a, const T & b)
{
#ifdef DEBUG_PSEUDO
    cout << "\nIN PSEUDODIVREM_reduce:  a = " << a << ", b = " << b << endl;
    T A = a;
    T B = b;
#endif

  register long i,j;
  long da, db, dq, LCIsOne;
  T lb;
  const S *bp;
  S *qp,*rp;
  S LC,temp;


  da = deg(a);
  db = deg(b);

  if (db < 0) Error("T: division by zero");

  if (da < db) {
    r = a;
    clear(q);
    set(d);
    return;
  }

  if (&q == &b) {
    lb = b;
    bp = lb.rep.elts();
  }
  else
    bp = b.rep.elts();

  if (IsOne(bp[db]))
    LCIsOne = 1;
  else {
    LCIsOne = 0;
    LC = bp[db];
  }


  // r is the remainder
  r = a;
  rp = r.rep.elts();
 
  // q is the quotient
  dq = da - db;
  q.rep.SetLength(dq+1);
  qp = q.rep.elts();


  set(d);

#ifdef DEBUG_PSEUDO
    cout << "LC = " << LC << endl;
    cout << "q = " << q << ", r = " << r << ", dr = " << dq+db << ", db = " << db << endl;
#endif

  for (j=dq; j >= 0; j--) {
    S t = rp[j+db];

    if (IsZero(t))
      clear(qp[j]);
    else {
      //
      // q = l(b) q + l(r) X^j, j = dr-db, 
      //
      if (!LCIsOne) {
	for (i=j+1; i<=dq; ++i)
	  mul(qp[i],qp[i],LC);
      }
      qp[j] = t;


      //
      // r = l(b) r - l(r) b X^j, j = dr-db
      //
      if (!LCIsOne) {
	for (i=0; i<j+db; ++i)
	  mul(rp[i],rp[i],LC);
      }

      for (i=0; i<db; ++i) {
	mul(temp,bp[i],t);
	sub(rp[i+j],rp[i+j],temp);
      }

      //
      // update multiplier
      //
      if (!LCIsOne)
        mul(d,d,LC);
    }

#ifdef DEBUG_PSEUDO
    cout << "j=" << j << ", q = " << q << ", r = " << r << ", dr = " << j+db << ", db = " << db << endl;
#endif
  }

  //
  // remove leading zeros of r
  //
  r.rep.SetLength(db);
  r.normalize();

#ifdef DEBUG_PSEUDO
    cout << "DONE:  d = " << d << ", q = " << q << ", r = " << r << endl;
    cout << "TEST:  d a = " << d*A << ", b q + r = " << B*q + r << endl;

    if (d*A != B*q + r) {
    cerr << "ERROR PSEUDO!" << endl;
    cerr << "a = " << A << ", b = " << B << endl;
    cerr << "d = " << d << ", q = " << q << ", r = " << r << endl;
    cerr << "d a = " << d*A << ", b q + r = " << B*q + r << endl;
    exit(1);
    }
#endif
}



template < class S, class T >
void PseudoDivRem_accumulate(S & d, T & q, T & r, const T & a, const T & b)
{
  // Calls PseudoDivRem_reduce by default.  Specializations exist for
  // ZZ_pX and GF2EX.

  PseudoDivRem_reduce(d,q,r,a,b);
}
