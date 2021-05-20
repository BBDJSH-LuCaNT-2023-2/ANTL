/**
 * @file xgcd_iter.hpp
 * @author Michael Jacobson
 *
 * @brief Routines for extended Euclidean algorithm but computing only one
 * of the Bezout coefficients (XGCD_LEFT) and partial extended Euclidean algorithm
 * that terminates when the remainder gets below a given threshold (XGCD_PARTIAL and
 * XGCD_PARTIAL_REDUCE). The main functions are:
 *   - XGCD_LEFT_ITER: extended Euclidean algorithm computing only one of the
 *     two solutions.  Specifically, compute G and X, given A, and B, such that
 *     A X = G (mod B).

 *   - XGCD_PARTIAL_ITER: extended Euclidean algorithm terminated when the size of
 *     a remainder goes below the given threshold.
 *
 *   - XGCD_PARTIAL_REDUCE_ITER: same as XGCD_PARTIAL_ITER, with an extra parameter
 *     for revised termination condition for ideal reduction algorithm
 *     (polynomial types only)
 *
 * These functions use an iterative method, either XGCD_*_NTL or XGCD_*_PSEUDO,
 * depending on the thresholds in qir_tune.cpp.
 *
 * Only defined for polynomial base types.
 */

#ifndef XGCD_ITER_H
#define XGCD_ITER_H

#include <ANTL/utilities.hpp>
#include <ANTL/thresholds.hpp>
#include <ANTL/XGCD/xgcd_plain.hpp>
#include <ANTL/XGCD/xgcd_pseudo.hpp>



//
// XGCD
//

template < class T >
void XGCD_ITER(T & G, T & X, T & Y, const T & A, const T & B);



//
// XGCD_LEFT
//

template < class T >
void XGCD_LEFT_ITER(T & G, T & X, const T & A, const T & B);



//
// Partial Euclidean algorithm (for NUCOMP)
//

template < class T >
void XGCD_PARTIAL_ITER(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag);

// flag not requred for GF2EX version
void XGCD_PARTIAL_ITER(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound);

// integer version uses a ZZ bound and Lehmer's variation (also no flag)
//void XGCD_PARTIAL_ITER(ZZ & R2, ZZ & R1, ZZ & C2, ZZ & C1, const ZZ & bound);



//
// Parital extended Euclidean algorithm (for fast reduce)
//

template < class T >
void XGCD_PARTIAL_REDUCE_ITER(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even);

// flag not requred for GF2EX version
void XGCD_PARTIAL_REDUCE_ITER(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even);



// Unspecialized template definitions.
#include "../../../src/XGCD/xgcd_iter_impl.hpp"

#endif // guard
