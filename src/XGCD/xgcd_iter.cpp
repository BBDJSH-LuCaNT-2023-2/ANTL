/**
 * @file xgcd_iter.cpp
 * @author Michael Jacobson
 * @remark specialized implementations of XGCD, XGCD_LEFT_ITER, XGCD_PARTIAL_ITER, and
 * XGCD_PARTIAL_REDUCE_ITER
 */

#include <ANTL/XGCD/xgcd_iter.hpp>

//
// Partial Euclidean algorithm (for NUCOMP)
//

void XGCD_PARTIAL_ITER(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound)
{
  if (deg(R1) < Thresholds<GF2EX>::get_pseudo_xgcd_partial_crossover(deg(R1),bound))
    XGCD_PARTIAL_PSEUDO(R2,R1,C2,C1,bound);
  else
    XGCD_PARTIAL_PLAIN(R2,R1,C2,C1,bound);
}



//
// Partial Euclidean algorithm (for fast reduce) with revised termination for
// polynomial types
//

void XGCD_PARTIAL_REDUCE_ITER(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even)
{
  if (deg(R1) < Thresholds<GF2EX>::get_pseudo_xgcd_partial_crossover(deg(R1),bound))
    XGCD_PARTIAL_REDUCE_PSEUDO(R2,R1,B2,B1,bound,even);
  else
    XGCD_PARTIAL_REDUCE_PLAIN(R2,R1,B2,B1,bound,even);
}
