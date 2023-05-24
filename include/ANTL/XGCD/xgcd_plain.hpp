/**
 * @file xgcd_plain.hpp
 * @author Michael Jacobson
 *
 * @brief Routines for extended Euclidean algorithm, extended Euclidean algorithm
 * computing only one of the Bezout coefficients (XGCD_LEFT), and partial extended
 * Euclidean algorithm that terminates when the remainder gets below a given threshold
 * (XGCD_PARTIAL and XGCD_PARTIAL_REDUCE). The main functions are:
 *   - XGCD_PLAIN: extended Euclidean algorithm computing only one of the
 *     two solutions.  Specifically, compute G, X, Y given A, and B, such that
 *     A X + B Y = G = gcd(A,B)
 *
 *   - XGCD_LEFT_PLAIN: extended Euclidean algorithm computing only one of the
 *     two solutions.  Specifically, compute G and X, given A, and B, such that
 *     A X = G (mod B).
 *
 *   - XGCD_PARTIAL_PLAIN: extended Euclidean algorithm terminated when the size of
 *     a remainder goes below the given threshold.  Uses partial extended Euclid for
 *     polynomials.
 *
 *   - XGCD_PARTIAL_REDUCE_PLAIN: same as XGCD_PARTIAL_NTL, with an extra parameter for
 *     revised termination condition for ideal reduction algorithm
 *     (polynomial types only)
 *
 * These functions use a basic version of the extended Eulidean algorithm, but modified when
 * possible (i.e., whenever the half gcd method is not used) to compute only one of
 * the two coefficients.
 *
 * Only defined for polynomial base types.
 */

#ifndef XGCD_PLAIN_H
#define XGCD_PLAIN_H


#include <ANTL/common.hpp>

// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.
NTL_CLIENT


//
// XGCD
//

template < class T >
void XGCD_PLAIN(T & G, T & X, T & Y, const T & A, const T & B);
template<>
void XGCD_PLAIN(long&, long&, long&, const long&, const long&);
template<>
void XGCD_PLAIN(ZZ&, ZZ&, ZZ&, const ZZ&, const ZZ&);



//
// XGCD_LEFT
//

template < class T >
void XGCD_LEFT_PLAIN(T & G, T & X, const T & A, const T & B);
template<>
void XGCD_LEFT_PLAIN(long & G, long & X, const long & A, const long & B);


//
// Partial Euclidean algorithm (for NUCOMP)
//

template < class T >
void XGCD_PARTIAL_PLAIN(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag);

// flag not requred for GF2EX version
void XGCD_PARTIAL_PLAIN(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound);

// integer version uses a ZZ bound and Lehmer's variation (also no flag)
//void XGCD_PARTIAL_PLAIN(ZZ & R2, ZZ & R1, ZZ & C2, ZZ & C1, const ZZ & bound);



//
// Parital extended Euclidean algorithm (for fast reduce)
//

template < class T >
void XGCD_PARTIAL_REDUCE_PLAIN(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even);

// flag not requred for GF2EX version
void XGCD_PARTIAL_REDUCE_PLAIN(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even);




// Unspecialized template definitions.
#include "../../../src/XGCD/xgcd_plain_impl.hpp"

#endif // guard
