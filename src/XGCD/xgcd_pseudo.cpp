/**
 * @file xgcd_pseudo.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of XGCD variations with pseudo-division
 */

#include <ANTL/XGCD/xgcd_pseudo.hpp>

//
// XGCD with one inversion
//

template <>
void XGCD_PSEUDO(GF2EX & G, GF2EX & X, GF2EX & Y, const GF2EX & A, const GF2EX & B)
{
  XGCD_PSEUDO_work<GF2E,GF2EX>(G,X,Y,A,B);
}

template <>
void XGCD_PSEUDO(zz_pX & G, zz_pX & X, zz_pX & Y, const zz_pX & A, const zz_pX & B)
{
  XGCD_PSEUDO_work<zz_p,zz_pX>(G,X,Y,A,B);
}

template <>
void XGCD_PSEUDO(ZZ_pX & G, ZZ_pX & X, ZZ_pX & Y, const ZZ_pX & A, const ZZ_pX & B)
{
  XGCD_PSEUDO_work<ZZ_p,ZZ_pX>(G,X,Y,A,B);
}

template <>
void XGCD_PSEUDO(zz_pEX & G, zz_pEX & X, zz_pEX & Y, const zz_pEX & A, const zz_pEX & B)
{
  XGCD_PSEUDO_work<zz_pE,zz_pEX>(G,X,Y,A,B);
}

template <>
void XGCD_PSEUDO(ZZ_pEX & G, ZZ_pEX & X, ZZ_pEX & Y, const ZZ_pEX & A, const ZZ_pEX & B)
{
  XGCD_PSEUDO_work<ZZ_pE,ZZ_pEX>(G,X,Y,A,B);
}



//
// XGCD_LEFT with one inversion
//

template <>
void XGCD_LEFT_PSEUDO(GF2EX & G, GF2EX & X, const GF2EX & A, const GF2EX & B)
{
  XGCD_LEFT_PSEUDO_work<GF2E,GF2EX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PSEUDO(zz_pX & G, zz_pX & X, const zz_pX & A, const zz_pX & B)
{
  XGCD_LEFT_PSEUDO_work<zz_p,zz_pX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PSEUDO(ZZ_pX & G, ZZ_pX & X, const ZZ_pX & A, const ZZ_pX & B)
{
  XGCD_LEFT_PSEUDO_work<ZZ_p,ZZ_pX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PSEUDO(zz_pEX & G, zz_pEX & X, const zz_pEX & A, const zz_pEX & B)
{
  XGCD_LEFT_PSEUDO_work<zz_pE,zz_pEX>(G,X,A,B);
}

template <>
void XGCD_LEFT_PSEUDO(ZZ_pEX & G, ZZ_pEX & X, const ZZ_pEX & A, const ZZ_pEX & B)
{
  XGCD_LEFT_PSEUDO_work<ZZ_pE,ZZ_pEX>(G,X,A,B);
}




//
// Partial Euclidean algorithm (for NUCOMP) with pseudo-division
//  - Assumes that R1 is reduced modulo R2
//

template <>
void XGCD_PARTIAL_PSEUDO(zz_pX & R2, zz_pX & R1, zz_pX & C2, zz_pX & C1, long bound, bool & flag)
{
  XGCD_PARTIAL_PSEUDO_work<zz_p,zz_pX>(R2,R1,C2,C1,bound,flag);
}

template <>
void XGCD_PARTIAL_PSEUDO(ZZ_pX & R2, ZZ_pX & R1, ZZ_pX & C2, ZZ_pX & C1, long bound, bool & flag)
{
  XGCD_PARTIAL_PSEUDO_work<ZZ_p,ZZ_pX>(R2,R1,C2,C1,bound,flag);
}

template <>
void XGCD_PARTIAL_PSEUDO(zz_pEX & R2, zz_pEX & R1, zz_pEX & C2, zz_pEX & C1, long bound, bool & flag)
{
  XGCD_PARTIAL_PSEUDO_work<zz_pE,zz_pEX>(R2,R1,C2,C1,bound,flag);
}

template <>
void XGCD_PARTIAL_PSEUDO(ZZ_pEX & R2, ZZ_pEX & R1, ZZ_pEX & C2, ZZ_pEX & C1, long bound, bool & flag)
{
  XGCD_PARTIAL_PSEUDO_work<ZZ_pE,ZZ_pEX>(R2,R1,C2,C1,bound,flag);
}



void 
XGCD_PARTIAL_PSEUDO(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound)
{
  GF2E d,z;

#ifdef DEBUG_PSEUDO
  cout << "\nIn XGCD_PARTIAL_PSEUDO(GF2EX)" << endl;
  cout << "R2 = " << R2 << endl;
  cout << "R1 = " << R1 << endl;
  GF2EX OR2 = R2;
#endif

  long e = deg(R2) + 1;
  GF2EX temp(INIT_SIZE, e), q(INIT_SIZE, e), r(INIT_SIZE,e),
    CC1(INIT_SIZE, e), CC2(INIT_SIZE, e);

  clear(CC2);
  set(CC1);

  set(z);

  while (deg (R1) > bound)
    {
      PseudoDivRem (d, q, R2, R2, R1);
      swap(R2,R1);
      mul(z,z,d);

      // r = d C2 + q C1
      mul(r,CC2,d);
      mul(temp,CC1,q);
      add(r,r,temp);
      CC2 = CC1;
      CC1 = r;
    }

#ifdef DEBUG_PSEUDO
  cout << "Done:  z = " << z << endl;
  cout << "R2 = " << R2 << endl;
  cout << "R1 = " << R1 << endl;
  cout << "C2 = " << CC2 << endl;
  cout << "C1 = " << CC1 << endl;
  cout << "C1 R2 + C2 R1 = " << C1*R2 + CC2*R1 << endl;
  cout << "OR2           = " << OR2 << endl;
  cout << "(C1 R2 + C2 R1) / OR2 = " << (C1*R2 + CC2*R1) / OR2 << endl;
#endif

  // remove inv(z)
  inv(z, z);
  mul(R1,R1,z);
  C2 = CC2;
  mul(C1,CC1,z);

#ifdef DEBUG_PSEUDO
  cout << "Norm:  z = " << z << endl;
  cout << "R2 = " << R2 << endl;
  cout << "R1 = " << R1 << endl;
  cout << "C2 = " << C2 << endl;
  cout << "C1 = " << C1 << endl;
  cout << "C1 R2 + C2 R1 = " << C1*R2 + CC2*R1 << endl;
  cout << "(C1 R2 + C2 R1) / OR2 = " << (C1*R2 + CC2*R1) / OR2 << endl;
  cout << "z = " << z << endl;
#endif
}



//
// Partial Euclidean algorithm (for fast reduce) with pseudo-division
//  - Assumes that R1 is reduced modulo R2
//

template <>
void XGCD_PARTIAL_REDUCE_PSEUDO(zz_pX & R2, zz_pX & R1, zz_pX & B2, zz_pX & B1, long bound, bool & flag, bool even)
{
  XGCD_PARTIAL_REDUCE_PSEUDO_work<zz_p,zz_pX>(R2,R1,B2,B1,bound,flag,even);
}

template <>
void XGCD_PARTIAL_REDUCE_PSEUDO(ZZ_pX & R2, ZZ_pX & R1, ZZ_pX & B2, ZZ_pX & B1, long bound, bool & flag, bool even)
{
  XGCD_PARTIAL_REDUCE_PSEUDO_work<ZZ_p,ZZ_pX>(R2,R1,B2,B1,bound,flag,even);
}

template <>
void XGCD_PARTIAL_REDUCE_PSEUDO(zz_pEX & R2, zz_pEX & R1, zz_pEX & B2, zz_pEX & B1, long bound, bool & flag, bool even)
{
  XGCD_PARTIAL_REDUCE_PSEUDO_work<zz_pE,zz_pEX>(R2,R1,B2,B1,bound,flag,even);
}

template <>
void XGCD_PARTIAL_REDUCE_PSEUDO(ZZ_pEX & R2, ZZ_pEX & R1, ZZ_pEX & B2, ZZ_pEX & B1, long bound, bool & flag, bool even)
{
  XGCD_PARTIAL_REDUCE_PSEUDO_work<ZZ_pE,ZZ_pEX>(R2,R1,B2,B1,bound,flag,even);
}




void 
XGCD_PARTIAL_REDUCE_PSEUDO(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even)
{
  GF2E d,z;

  long e = deg(R2) + 1;
  GF2EX temp(INIT_SIZE, e), q(INIT_SIZE, e), r(INIT_SIZE,e),
    BB1(INIT_SIZE, e), BB2(INIT_SIZE, e);

  set(BB2);
  clear(BB1);

  set(z);

  while (deg (R1) > bound || (deg(R1) == bound && !even))
     {
      PseudoDivRem (d, q, R2, R2, R1);
      swap(R2,R1);
      mul(z,z,d);

      // r = d B2 + q B1
      mul(r,BB2,d);
      mul(temp,BB1,q);
      add(r,r,temp);
      BB2 = BB1;
      BB1 = r;
    }

  // remove inv(z)
  inv(z, z);
  mul(R1,R1,z);
  B2 = BB2;
  mul(B1,BB1,z);
}
