/**
 * @file hxgcd.cpp
 * @author Laurent Imbert
 * @remark
 */


// NOT SURE THESE FUNCTIONS ARE EVER CALLED!!!!

#include <ANTL/XGCD/hxgcd.hpp>
#include <ANTL/XGCD/xgcd_iter.hpp>

template <>
void HXGCD(ZZ_pX& G, ZZ_pX& U, ZZ_pX& V, const ZZ_pX& A, const ZZ_pX& B)
{
  // Uses NTL XGCD implementation
  XGCD(G, U, V, A, B);
}


// -----------------------------------------------------------------------------

template <>
void HXGCD(zz_pX& G, zz_pX& U, zz_pX& V, const zz_pX& A, const zz_pX& B)
{
  // Uses NTL XGCD implementation
  XGCD(G, U, V, A, B);
}


//
// Partial Euclidean algorithm (for NUCOMP and fast reduce)
//

void HXGCD_PARTIAL(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound)
{
  TMatrix<GF2EX> M_out;

  // if deg(R2) < deg(R1), first quotient is zero.  Swap.
  bool swapped = false;
  if (deg(R2) < deg(R1)) {
    swap(R2,R1);
    swapped = true;
  }

  long d_red = deg(R2) - bound;
  HALF_GCD(M_out, R2, R1, d_red);

  // swap B2 and B1, because of initial swap
  if (swapped) {
    C2 = M_out(1,1);
    C1 = M_out(0,1);
  }
  else {
    C2 = M_out(0,1);
    C1 = M_out(1,1);
  }
}



//
// Partial Euclidean algorithm (for fast reduce) with revised termination for
// polynomial types
//

void HXGCD_PARTIAL_REDUCE(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even)
{
  TMatrix<GF2EX> M_out;

  // if deg(R2) < deg(R1), first quotient is zero.  Swap.
  bool swapped = false;
  if (deg(R2) < deg(R1)) {
    swap(R2,R1);
    swapped = true;
  }

  long d_red = deg(R2) - bound;

  HALF_GCD(M_out, R2, R1, d_red);

  // swap B2 and B1, because of initial swap
  if (swapped) {
    B2 = M_out(0,1);
    B1 = M_out(1,1);
  }
  else {
    B2 = M_out(0,0);
    B1 = M_out(1,0);
  }

  if (deg(R1) == bound && !even) {
    GF2EX Q,t;

    DivRem(Q, R2, R2, R1);
    swap(R1, R2);

      // r = B2 + q B1
      mul(t,Q,B1);
      add(t,B2,t);
      B2 = B1;
      B1 = t;
  }
}
