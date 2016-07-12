/**
 * @file pseudodiv.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of pseudodivision.
 */

#include <ANTL/Arithmetic/pseudodiv.hpp>

//
// Pseudo division with remainder (for 1 inversion XGCD routines)
//

// Pseudodivision:  zz_pX version

template <>
void PseudoDivRem(zz_p & d, zz_pX & q, zz_pX & r, const zz_pX & a, const zz_pX & b)
{
#ifdef DEBUG_PSEUDO
  cout << "\nIN PSEUDODIVREM(zz_pX):  a = " << a << ", b = " << b << endl;
  zz_pX A = a;
  zz_pX B = b;
#endif

  register long i,j;
  long da, db, dq, LCIsOne;
  zz_pX lb;
  const zz_p *bp;
  zz_p *qp, *xp;

  zz_p t;
  long LCrep=0,S,T;
  mulmod_precon_t LCpinv=0;

  da = deg(a);
  db = deg(b);

  if (db < 0) Error("PSEUDODIVREM(zz_pX): division by zero");

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
    LCrep = rep(bp[db]);
  }


  // x is the remainder
  vec_zz_p x;
  if (&r == &a)
    xp = r.rep.elts();
  else {
    x = a.rep;
    xp = x.elts();
  }


  // q is the quotient
  dq = da - db;
  q.rep.SetLength(dq+1);
  qp = q.rep.elts();


  // set up modulus for mod p arithmetic
  long p = zz_p::modulus();
  mulmod_t pinv = zz_p::ModulusInverse();

  // set up preconditioner for modular multiplication by LC
  if (!LCIsOne)
    LCpinv = PrepMulModPrecon(LCrep, p, pinv); // ((double) T)*pinv;

#ifdef DEBUG_PSEUDO
  cout << "LC = " << LCrep << endl;
  cout << "q = " << q << ", r = " << r << ", dr = " << dq+db << ", db = " << db << endl;
#endif

  set(d);

  for (j=dq; j >= 0; j--) {
    t = xp[j+db];

    if (IsZero(t))
      clear(qp[j]);
    else {
      //
      // q = l(b) q + l(r) X^j, j = dr-db, 
      //
      if (!LCIsOne) {
	for (i=j+1; i<=dq; ++i)
	  qp[i].LoopHole() = MulModPrecon(rep(qp[i]),LCrep,p,LCpinv);
      }
      qp[j] = t;


      //
      // r = l(b) r - l(r) b X^j, j = dr-db
      //
      if (!LCIsOne) {
	for (i=0; i<j+db; ++i)
	  xp[i].LoopHole() = MulModPrecon(rep(xp[i]),LCrep,p,LCpinv);
      }

      T = rep(t);
      mulmod_precon_t Tpinv = PrepMulModPrecon(T, p, pinv); // ((double) T)*pinv;

      for (i=0; i<db; ++i) {
	S = MulModPrecon(rep(bp[i]),T,p,Tpinv);
	xp[i+j].LoopHole() = SubMod(rep(xp[i+j]),S,p);
      }

      //
      // update multiplier
      //
      if (!LCIsOne)
	d.LoopHole() = MulModPrecon(rep(d),LCrep,p,LCpinv);
    }

#ifdef DEBUG_PSEUDO
    cout << "j=" << j << ", q = " << q << ", r = " << r << ", dr = " << j+db << ", db = " << db << endl;
#endif
  }


  //
  // remove leading zeros of r
  //
  r.rep.SetLength(db);
  if (&r != &a) {
    for (i = 0; i < db; i++)
      r.rep[i] = xp[i];
  }
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



template <>
void PseudoDivRem_reduce(zz_p & d, zz_pX & q, zz_pX & r, const zz_pX & a, const zz_pX & b)
{
  PseudoDivRem(d,q,r,a,b);
}

template <>
void PseudoDivRem_accumulate(zz_p & d, zz_pX & q, zz_pX & r, const zz_pX & a, const zz_pX & b)
{
  PseudoDivRem(d,q,r,a,b);
}




//
// Pseudodivision:  ZZ_pX versions
//

template <>
void PseudoDivRem(ZZ_p & d, ZZ_pX & q, ZZ_pX & r, const ZZ_pX & a, const ZZ_pX & b)
{
  /*
  if (deg(b) > QIR_ZZ_pX_PSEUDODIV_CROSSOVER && deg(a) - deg(b) > QIR_ZZ_pX_PSEUDODIV_CROSSOVER)
    PseudoDivRem_reduce(d,q,r,a,b);
  else
    PseudoDivRem_accumulate(d,q,r,a,b);
  */
  PseudoDivRem_accumulate(d,q,r,a,b);
}



template <>
void PseudoDivRem_accumulate(ZZ_p & d, ZZ_pX & q, ZZ_pX & r, const ZZ_pX & a, const ZZ_pX & b)
{
#ifdef DEBUG_PSEUDO
  cout << "\nIN PSEUDODIVREM(ZZ_pX):  a = " << a << ", b = " << b << endl;
  ZZ_pX A = a;
  ZZ_pX B = b;
#endif

  register long i,j;
  long da, db, dq, LCIsOne;
  ZZ_pX lb;
  const ZZ_p *bp;
  ZZ_p *qp;

  ZZ_p LC,t;
  ZZ LCrep,s, *xp;

  da = deg(a);
  db = deg(b);

  if (db < 0) Error("PSEUDODIVREM(ZZ_pX): division by zero");

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
    LCrep = rep(bp[db]);
  }


  // q is the quotient
  dq = da - db;
  q.rep.SetLength(dq+1);
  qp = q.rep.elts();


  // x is the remainder
  ZZVec x(da + 1, ZZ_pInfo->ExtendedModulusSize + (dq+1)*ZZ_pInfo->size);
  for (i=0; i<=da; i++)
    x[i] = rep(a.rep[i]);
  xp = x.elts();


  set(d);

#ifdef DEBUG_PSEUDO
  cout << "LC = " << LC << endl;
  cout << "q = " << q << ", r = " << r << ", dr = " << dq+db << ", db = " << db << endl;
#endif

  for (j=dq; j >= 0; j--) {
    conv(t,xp[j+db]);

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
	  mul(xp[i],xp[i],LCrep);
      }

      for (i=0; i<db; ++i)
	MulSubFrom(xp[i+j],rep(bp[i]),rep(t));

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
  for (i = 0; i < db; i++)
    conv(r.rep[i],xp[i]);
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



//
// Pseudodivision:  GF2EX version
//

template <>
void PseudoDivRem(GF2E & d, GF2EX & q, GF2EX & r, const GF2EX & a, const GF2EX & b)
{
  /*
  if (deg(b) > QIR_GF2EX_PSEUDODIV_CROSSOVER && deg(a) - deg(b) > QIR_GF2EX_PSEUDODIV_CROSSOVER)
    PseudoDivRem_reduce(d,q,r,a,b);
  else
    PseudoDivRem_accumulate(d,q,r,a,b);
  */
//  PseudoDivRem_accumulate(d,q,r,a,b);
    PseudoDivRem_reduce(d,q,r,a,b);
}



template <>
void PseudoDivRem_accumulate(GF2E & d, GF2EX & q, GF2EX & r, const GF2EX & a, const GF2EX & b)
{
#ifdef DEBUG_PSEUDO
  cout << "\nIN PSEUDODIVREM(GF2EX):  a = " << a << ", b = " << b << endl;
  GF2EX A = a;
  GF2EX B = b;
#endif

  register long i,j;
  long da, db, dq, LCIsOne;
  GF2EX lb;
  const GF2E *bp;
  GF2E *qp;

  GF2E LC,t;
  GF2X LCrep,s, *xp;


  da = deg(a);
  db = deg(b);

  if (db < 0) Error("PSEUDODIVREM(GF2EX): division by zero");

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
    LCrep = rep(bp[db]);
  }


  // q is the quotient
  dq = da - db;
  q.rep.SetLength(dq+1);
  qp = q.rep.elts();


  // x is the remainder
  GF2XVec x(da + 1, (dq+2)*GF2E::WordLength());
  for (i=0; i<=da; i++)
    x[i] = rep(a.rep[i]);
  xp = x.elts();


  set(d);

#ifdef DEBUG_PSEUDO
  cout << "LC = " << LC << endl;
  cout << "q = " << q << ", r = " << r << ", dr = " << dq+db << ", db = " << db << endl;
#endif

  for (j=dq; j >= 0; j--) {
    conv(t,xp[j+db]);

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
  	  mul(xp[i],xp[i],LCrep);
      }

      for (i=0; i<db; ++i) {
	mul(s,rep(bp[i]),rep(t));
	add(xp[i+j],xp[i+j],s);
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
  for (i = 0; i < db; i++)
    conv(r.rep[i],xp[i]);
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
