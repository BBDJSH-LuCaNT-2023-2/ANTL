/**
 * @file mul_exact.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of multiply functions computing only high-degree
 * terms.
 */

#include <ANTL/Arithmetic/mul_exact.hpp>

//
// PlainMulExact
//

template <>
void
PlainMulExact(GF2EX & x, const GF2EX & a, const GF2EX & b, long n)
{
   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   if (&a == &b) {
     SqrExact(x, a, n);
      return;
   }

   long d = da+db;

   const GF2E *ap, *bp;
   GF2E *xp;
   
   GF2EX la, lb;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();

   if (&x == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   GF2X t, accum;

   for (i = n; i <= d; i++) {
      jmin = max(0, i-db);
      jmax = min(da, i);
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(bp[i-j]));
	 add(accum, accum, t);
      }
      conv(xp[i], accum);
   }
   x.normalize();
}


template <>
void
PlainMulExact(zz_pX& x, const zz_pX& a, const zz_pX& b, long n)
{
   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   if (da == 0) {
      mul(x, b, a.rep[0]);
      return;
   }

   if (db == 0) {
      mul(x, a, b.rep[0]);
      return;
   }

   long d = da+db;



   const zz_p *ap, *bp;
   zz_p *xp;
   
   zz_pX la, lb;

   if (&x == &a) {
     la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();

   if (&x == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   x.rep.SetLength(d+1);
   xp = x.rep.elts();

   long i, j, jmin, jmax;
   //   static long accum;
   zz_p t;

   for (i = n; i <= d; i++) {
      jmin = max(0, i-db);
      jmax = min(da, i);
      clear(xp[i]);
      for (j = jmin; j <= jmax; j++) {
	mul(t,ap[j],bp[i-j]);
	add(xp[i],xp[i],t);
      }
      //      conv(xp[i], accum);
   }
   x.normalize();
}


template <>
void
PlainMulExact(ZZ_pX& x, const ZZ_pX& a, const ZZ_pX& b, long n)
{
   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   if (da == 0) {
      mul(x, b, a.rep[0]);
      return;
   }

   if (db == 0) {
      mul(x, a, b.rep[0]);
      return;
   }

   long d = da+db;



   const ZZ_p *ap, *bp;
   ZZ_p *xp;
   
   ZZ_pX la, lb;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();

   if (&x == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   static ZZ t, accum;

   for (i = n; i <= d; i++) {
      jmin = max(0, i-db);
      jmax = min(da, i);
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(bp[i-j]));
	 add(accum, accum, t);
      }
      conv(xp[i], accum);
   }
   x.normalize();
}


template <>
void
PlainMulExact(ZZ_pEX& x, const ZZ_pEX& a, const ZZ_pEX& b, long n)
{
   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   long d = da+db;



   const ZZ_pE *ap, *bp;
   ZZ_pE *xp;
   
   ZZ_pEX la, lb;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();

   if (&x == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   ZZ_pX t, accum;

   for (i = n; i <= d; i++) {
      jmin = max(0, i-db);
      jmax = min(da, i);
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(bp[i-j]));
	 add(accum, accum, t);
      }
      conv(xp[i], accum);
   }
   x.normalize();
}


template <>
void
PlainMulExact(zz_pEX& x, const zz_pEX& a, const zz_pEX& b, long n)
{
   long da = deg(a);
   long db = deg(b);

   if (da < 0 || db < 0) {
      clear(x);
      return;
   }

   long d = da+db;



   const zz_pE *ap, *bp;
   zz_pE *xp;
   
   zz_pEX la, lb;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();

   if (&x == &b) {
      lb = b;
      bp = lb.rep.elts();
   }
   else
      bp = b.rep.elts();

   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   zz_pX t, accum;

   for (i = n; i <= d; i++) {
      jmin = max(0, i-db);
      jmax = min(da, i);
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(bp[i-j]));
	 add(accum, accum, t);
      }
      conv(xp[i], accum);
   }
   x.normalize();
}




//
// PlainSqrExact
//

template <>
void
PlainSqrExact(GF2EX & x, const GF2EX & a, long n)
{
   long da = deg(a);

   if (da < 0) {
      clear(x);
      return;
   }

   x.rep.SetLength(2*da+1);
   long i;

   long bnd = (n+1) >> 1;
   for (i = da; i >= bnd; i--) {
      sqr(x.rep[2*i], a.rep[i]);
      clear(x.rep[2*i-1]);
   }

   x.normalize();
}


template <>
void
PlainSqrExact(ZZ_pX& x, const ZZ_pX& a, long n)
{
   long da = deg(a);

   if (da < 0) {
      clear(x);
      return;
   }

   long d = 2*da;

   const ZZ_p *ap;
   ZZ_p *xp;

   ZZ_pX la;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();


   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   long m, m2;
   static ZZ t, accum;

   for (i = n; i <= d; i++) {
      jmin = max(0, i-da);
      jmax = min(da, i);
      m = jmax - jmin + 1;
      m2 = m >> 1;
      jmax = jmin + m2 - 1;
      clear(accum);
      for (j = jmin; j <= jmax; j++) {
	 mul(t, rep(ap[j]), rep(ap[i-j]));
	 add(accum, accum, t);
      }
      add(accum, accum, accum);
      if (m & 1) {
	 sqr(t, rep(ap[jmax + 1]));
	 add(accum, accum, t);
      }

      conv(xp[i], accum);
   }

   x.normalize();
}


template <>
void 
PlainSqrExact(zz_pX& x, const zz_pX& a, long n)
{
   long da = deg(a);

   if (da < 0) {
      clear(x);
      return;
   }

   long d = 2*da;

   const zz_p *ap;
   zz_p *xp;

   zz_pX la;

   if (&x == &a) {
      la = a;
      ap = la.rep.elts();
   }
   else
      ap = a.rep.elts();


   x.rep.SetLength(d+1);

   xp = x.rep.elts();

   long i, j, jmin, jmax;
   long m, m2;
   //   static long t, accum;
   zz_p t;
   
   for (i = n; i <= d; i++) {
      jmin = max(0, i-da);
      jmax = min(da, i);
      m = jmax - jmin + 1;
      m2 = m >> 1;
      jmax = jmin + m2 - 1;
      clear(xp[i]); 
      for (j = jmin; j <= jmax; j++) {
	mul(t,ap[j],ap[i-j]);
	add(xp[i],xp[i],t);
      }
      add(xp[i],xp[i],xp[i]);
      if (m & 1) {
         t = ap[jmax + 1];
	 sqr(t,ap[jmax+1]);
	 add(xp[i],xp[i],t);
      }

      //      conv(xp[i], accum);
   }

   x.normalize();
}


template <>
void
PlainSqrExact(ZZ_pEX& x, const ZZ_pEX& a, long n)
{
   if (IsZero(a)) {
      clear(x);
      return;
   }

   if (deg(a) == 0) {
      ZZ_pE res;
      sqr(res, ConstTerm(a));
      conv(x, res);
      return;
   } 

   // general case...Kronecker subst

   ZZ_pX A, C;

   long da = deg(a);

   long nn = ZZ_pE::degree();
   long n2 = 2*nn-1;

   if (NTL_OVERFLOW(2*da+1, n2, 0))
      Error("overflow in ZZ_pEX sqr");

   long i, j;

   A.rep.SetLength((da+1)*n2);

   for (i = 0; i <= da; i++) {
      const ZZ_pX& coeff = rep(a.rep[i]);
      long dcoeff = deg(coeff);
      for (j = 0; j <= dcoeff; j++)
         A.rep[n2*i + j] = coeff.rep[j]; 
   }

   A.normalize();

   sqr(C, A);

   long Clen = C.rep.length();
   long lc = (Clen + n2 - 1)/n2;
   long dc = lc - 1;

   x.rep.SetLength(dc+1);

   ZZ_pX tmp;
   
   for (i = n; i <= dc; i++) {
      tmp.rep.SetLength(n2);
      for (j = 0; j < n2 && n2*i + j < Clen; j++)
         tmp.rep[j] = C.rep[n2*i + j];
      for (; j < n2; j++)
         clear(tmp.rep[j]);
      tmp.normalize();
      conv(x.rep[i], tmp);
   }
  
  
   x.normalize();
}


template <>
void
PlainSqrExact(zz_pEX& x, const zz_pEX& a, long n)
{
   if (IsZero(a)) {
      clear(x);
      return;
   }

   if (deg(a) == 0) {
      zz_pE res;
      sqr(res, ConstTerm(a));
      conv(x, res);
      return;
   } 

   // general case...Kronecker subst

   zz_pX A, C;

   long da = deg(a);

   long nn = zz_pE::degree();
   long n2 = 2*nn-1;

   if (NTL_OVERFLOW(2*da+1, n2, 0))
      Error("overflow in zz_pEX sqr");

   long i, j;

   A.rep.SetLength((da+1)*n2);

   for (i = 0; i <= da; i++) {
      const zz_pX& coeff = rep(a.rep[i]);
      long dcoeff = deg(coeff);
      for (j = 0; j <= dcoeff; j++)
         A.rep[n2*i + j] = coeff.rep[j]; 
   }

   A.normalize();

   sqr(C, A);

   long Clen = C.rep.length();
   long lc = (Clen + n2 - 1)/n2;
   long dc = lc - 1;

   x.rep.SetLength(dc+1);

   zz_pX tmp;
   
   for (i = n; i <= dc; i++) {
      tmp.rep.SetLength(n2);
      for (j = 0; j < n2 && n2*i + j < Clen; j++)
         tmp.rep[j] = C.rep[n2*i + j];
      for (; j < n2; j++)
         clear(tmp.rep[j]);
      tmp.normalize();
      conv(x.rep[i], tmp);
   }
  
  
   x.normalize();
}


//
// SqrExact
//

/*
template <>
void
SqrExact(GF2EX & x, const GF2EX & a, long n)
{
   long da = deg(a);

   if (da < 0) {
      clear(x);
      return;
   }

   x.rep.SetLength(2*da+1);
   long i;

   long bnd = (n+1) >> 1;
   for (i = da; i >= bnd; i--) {
      sqr(x.rep[2*i], a.rep[i]);
      clear(x.rep[2*i-1]);
   }

   x.normalize();
}
*/


// Truncated Karatsuba multiplication

/*
static 
void PlainMul1(GF2X *xp, const GF2X *ap, long sa, const GF2X& b)
{
   long i;

   for (i = 0; i < sa; i++)
      mul(xp[i], ap[i], b);
}




static inline
void q_add(GF2X& x, const GF2X& a, const GF2X& b)

// This is a quick-and-dirty add rotine used by the karatsuba routine.
// It assumes that the output already has enough space allocated,
// thus avoiding any procedure calls.
// WARNING: it also accesses the underlying WordVector representation
// directly...that is dirty!.
// It shaves a few percent off the running time.

{
   _ntl_ulong *xp = x.xrep.elts();
   const _ntl_ulong *ap = a.xrep.elts();
   const _ntl_ulong *bp = b.xrep.elts();

   long sa = ap[-1];
   long sb = bp[-1];

   long i;

   if (sa == sb) {
      for (i = 0; i < sa; i++)
         xp[i] = ap[i] ^ bp[i];

      i = sa-1;
      while (i >= 0 && !xp[i]) i--;
      xp[-1] = i+1;
   }
   else if (sa < sb) {
      for (i = 0; i < sa; i++)
         xp[i] = ap[i] ^ bp[i];

      for (; i < sb; i++)
         xp[i] = bp[i];

      xp[-1] = sb;
   }
   else { // sa > sb
      for (i = 0; i < sb; i++)
         xp[i] = ap[i] ^ bp[i];

      for (; i < sa; i++)
         xp[i] = ap[i];

      xp[-1] = sa;
   }
}


static inline
void q_copy(GF2X& x, const GF2X& a)
// see comments for q_add above

{
   _ntl_ulong *xp = x.xrep.elts();
   const _ntl_ulong *ap = a.xrep.elts();

   long sa = ap[-1];
   long i;

   for (i = 0; i < sa; i++)
      xp[i] = ap[i];

   xp[-1] = sa;
}



static
void KarFold(GF2X *T, const GF2X *b, long sb, long hsa)
{
   long m = sb - hsa;
   long i;

   for (i = 0; i < m; i++)
      q_add(T[i], b[i], b[hsa+i]);

   for (i = m; i < hsa; i++)
      q_copy(T[i], b[i]);
}


static
void KarAdd(GF2X *T, const GF2X *b, long sb)
{
   long i;

   for (i = 0; i < sb; i++)
      q_add(T[i], T[i], b[i]);
}

static
void KarFix(GF2X *c, const GF2X *b, long sb, long hsa)
{
   long i;

   for (i = 0; i < hsa; i++)
      q_copy(c[i], b[i]);

   for (i = hsa; i < sb; i++)
      q_add(c[i], c[i], b[i]);
}




static
void KarMulExactRec(GF2X *c, const GF2X *a, 
            long sa, const GF2X *b, long sb, GF2X *stk, long offset, long n)
{
#ifdef DEBUG_MULEXACT
  //  cout << "REC:  sa = " << sa << ", sb = " << sb << ", offset = " << offset << ", n = " << n << endl;
#endif

   if (sb == 1) {  
      if (sa == 1) 
         mul(*c, *a, *b);
      else
	PlainMul1(c, a, sa, *b);

      return;
   }

   if (sb == 2 && sa == 2) {
      mul(c[0], a[0], b[0]);
      mul(c[2], a[1], b[1]);
      q_add(stk[0], a[0], a[1]);
      q_add(stk[1], b[0], b[1]);
      mul(c[1], stk[0], stk[1]);
      q_add(c[1], c[1], c[0]);
      q_add(c[1], c[1], c[2]);
      
      return;
   }

   long hsa = (sa + 1) >> 1;

#ifdef DEBUG_MULEXACT
   //   cout << "hsa = " << hsa << endl;
#endif

   if (hsa < sb) {
     // normal case

      long hsa2 = hsa << 1;
      long off2 = offset + hsa + sa;

#ifdef DEBUG_MULEXACT
      //      cout << "off2 = " << off2 << endl;
#endif

      // recursively compute a_hi * b_hi into high part of c
      KarMulExactRec(c + hsa2, a+hsa, sa-hsa, b+hsa, sb-hsa, stk, offset + hsa2, n);

      if (off2 >= n) {
	GF2X *T1, *T2, *T3;

	T1 = stk; stk += hsa;
	T2 = stk; stk += hsa;
	T3 = stk; stk += hsa2 - 1;

	// compute T1 = a_lo + a_hi

	KarFold(T1, a, sa, hsa);

	// compute T2 = b_lo + b_hi

	KarFold(T2, b, sb, hsa);

	// recursively compute T3 = T1 * T2

	KarMulExactRec(T3, T1, hsa, T2, hsa, stk, offset+hsa ,n);

        // subtract a_hi * b_hi from T3
	KarAdd(T3, c + hsa2, sa + sb - hsa2 - 1);


        // recursively compute a_lo*b_lo into low part of c
        // and subtract from T3
 	KarMulExactRec(c, a, hsa, b, hsa, stk, offset+hsa, n);
	KarAdd(T3, c, hsa2 - 1);

	clear(c[hsa2 - 1]);

	// finally, add T3 * X^{hsa} to c

	KarAdd(c+hsa, T3, hsa2-1);
      }
   }
   else {
     // degenerate case

      GF2X *T;

      T = stk; stk += hsa + sb - 1;

      // recursively compute b*a_hi into high part of c

      KarMulExactRec(c + hsa, a + hsa, sa - hsa, b, sb, stk, offset + hsa, n);

      // recursively compute b*a_lo into T
      long off2 = offset + hsa + sb;
      if (off2 >= n) {
	KarMulExactRec(T, a, hsa, b, sb, stk, offset, n);
	KarFix(c, T, hsa + sb - 1, hsa);
      }
   }
}




template <>
void  qir_mul_exact<GF2EX>::
KarMulExact(GF2EX & x, const GF2EX & a, const GF2EX & b, long n)
{
#ifdef DEBUG_MULEXACT
  cout << "\nKarMulExact" << endl;
  cout << "a = " << a << endl;
  cout << "b = " << b << endl;
  cout << "n = " << n << endl;
  GF2EX pp = a*b;
  GF2EX ans = a*b - NTL::MulExact(a,b,n);
  cout << "a*b = " << pp << endl;
  cout << "ans = " << ans << endl;
#endif

   long dn, hn, sp;
   GF2EX aa=a,bb=b;

   long sa = aa.rep.length();
   long sb = bb.rep.length();

   if (sa < sb) {
     GF2EX t=aa; aa=bb; bb=t;
     long tl=sa; sa=sb; sb=tl;
   }

   dn = max(sa, sb);
   sp = 0;
   do {
      hn = (dn+1) >> 1;
      sp += (hn << 2) - 1;
      dn = hn;
   } while (dn > 1);

   GF2XVec stk;
   stk.SetSize(sp + 2*(sa+sb)-1, 2*GF2E::WordLength()); 

   long i;

   for (i = 0; i < sa; i++)
      stk[i+sa+sb-1] = rep(aa.rep[i]);

   for (i = 0; i < sb; i++)
      stk[i+2*sa+sb-1] = rep(bb.rep[i]);

   KarMulExactRec(&stk[0], &stk[sa+sb-1], sa, &stk[2*sa+sb-1], sb, 
		  &stk[2*(sa+sb)-1],0,n);

   x.rep.SetLength(sa+sb-1);

   for (i = 0; i < sa+sb-1; i++)
      conv(x.rep[i], stk[i]);

   x.normalize();

#ifdef DEBUG_MULEXACT
   cout << "x = " << x << endl;
   GF2EX xt = x - trunc(x,n);
   cout << "xt = " << xt << endl;
   if (ans != xt) {
     cout << "ERROR KarMulExact" << endl;
     exit(1);
   }
#endif
}
*/
