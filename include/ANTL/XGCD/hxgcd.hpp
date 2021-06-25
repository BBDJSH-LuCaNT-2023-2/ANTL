/**
 * @file hxgcd.hpp
 * @author Laurent Imbert
 * @brief  Routines for extended Euclidean algorithm using the half-gcd method
 */

#ifndef HXGCD_H
#define HXGCD_H

#include <ANTL/common.hpp>
#include <ANTL/thresholds.hpp>
#include <ANTL/XGCD/xgcd_iter.hpp>


template < class T >
void HXGCD_LEFT(T & G, T & X, const T & A, const T & B);

template < class T >
void HXGCD(T & G, T & X, T & Y, const T & A, const T & B);



//
// Partial Euclidean algorithm (for NUCOMP)
//

template < class T >
void HXGCD_PARTIAL(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag);

// flag not requred for GF2EX version
void HXGCD_PARTIAL(GF2EX & R2, GF2EX & R1, GF2EX & C2, GF2EX & C1, long bound);



//
// Parital extended Euclidean algorithm (for fast reduce)
//

template < class T >
void HXGCD_PARTIAL_REDUCE(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even);

// flag not requred for GF2EX version
void HXGCD_PARTIAL_REDUCE(GF2EX & R2, GF2EX & R1, GF2EX & B2, GF2EX & B1, long bound, bool even);



//
// Declare specialized methods
//

template <> void HXGCD(ZZ_pX& G, ZZ_pX& U, ZZ_pX& V, const ZZ_pX& A, const ZZ_pX& B);
template <> void HXGCD(zz_pX& G, zz_pX& U, zz_pX& V, const zz_pX& A, const zz_pX& B);



// Unspecialized template definitions.
#include "../../../src/XGCD/hxgcd_impl.hpp"

#endif // guard
