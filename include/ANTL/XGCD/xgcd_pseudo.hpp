/**
 * @file xgcd_pseudo.hpp
 * @author Michael Jacobson
 *
 * @brief Routines for partial extended Euclidean algorithm, where the
 * algorithm terminates once the size of the remainder sequence falls below
 * an input threshold.  The main functions are:
 *   - XGCD_LEFT_PSEUDO:  uses pseudodivision so that only one field inversion
 *     is required (for polynomial operands only)
 *
 *   - XGCD_PARTIAL_PSEUDO: same as XGCD_PARTIAL, but uses pseudodivision so
 *     so that only one field inverse is required (poly types)
 * 
 *   - XGCD_PARTIAL_REDUCE_PSEUDO: same as XGCD_PARTIAL_PSEUDO, plus extra
 *     parameter for termination condition for reduction algorithm
 *
 * Only defined for polynomial base types.
 */

#ifndef XGCD_PSEUDO_H
#define XGCD_PSEUDO_H

#include <ANTL/utilities.hpp>
#include <ANTL/Arithmetic/pseudodiv.hpp>


// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.

NTL_CLIENT


//
// XGCD with only one field inversion

template < class T >
void XGCD_PSEUDO(T & G, T & X, T & Y, const T & A, const T & B);


//
// XGCD_LEFT with only one field inversion

template < class T >
void XGCD_LEFT_PSEUDO(T & G, T & X, const T & A, const T & B);


//
// Partial Euclidean algorithm (for NUCOMP) using pseudodivision
//

template < class T >
void XGCD_PARTIAL_PSEUDO(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag);

// flag not requred for GF2EX version
void XGCD_PARTIAL_PSEUDO(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound);


//
// Partial Euclidean algorithm (for fast reduce) using pseudodivision
//

template < class T >
void XGCD_PARTIAL_REDUCE_PSEUDO(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even);

// flag not requred for GF2EX version
void XGCD_PARTIAL_REDUCE_PSEUDO(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even);



// Unspecialized template definitions.
#include "../src/XGCD/xgcd_pseudo_impl.hpp"

#endif // guard
