/**
 * @file mul_exact.hpp
 * @author Michael Jacobson
 * @brief methods for multiplication in which only the upper degree terms
 * are computed (used for subsequent exact divisions)
 */

#ifndef MUL_EXACT_H
#define MUL_EXACT_H

#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ.h>
#include <ANTL/thresholds.hpp>

NTL_CLIENT


template < class T >
void PlainMulExact(T & x, const T & a, const T & b, long n);

template < class T >
void PlainSqrExact(T & x, const T & a, long n);

template < class T >
void MulExact(T & x, const T & a, const T & b, long n);

template < class T >
void SqrExact(T & x, const T & a, long n);


//
// Declare specialized methods
//

template <> void PlainMulExact(GF2EX & x, const GF2EX & a, const GF2EX & b, long n);
template <> void PlainMulExact(ZZ_pX & x, const ZZ_pX & a, const ZZ_pX & b, long n); 
template <> void PlainMulExact(zz_pX & x, const zz_pX & a, const zz_pX & b, long n);
template <> void PlainMulExact(ZZ_pEX & x, const ZZ_pEX & a, const ZZ_pEX & b, long n);
template <> void PlainMulExact(zz_pEX & x, const zz_pEX & a, const zz_pEX & b, long n);

template <> void PlainSqrExact(GF2EX & x, const GF2EX & a, long n);
template <> void PlainSqrExact(ZZ_pX & x, const ZZ_pX & a, long n);
template <> void PlainSqrExact(zz_pX & x, const zz_pX & a, long n);
template <> void PlainSqrExact(ZZ_pEX & x, const ZZ_pEX & a, long n);
template <> void PlainSqrExact(zz_pEX & x, const zz_pEX & a, long n);

//template <> void SqrExact(GF2EX & x, const GF2EX & a, long n);


// Unspecialized template definitions.
#include "../src/Arithmetic/mul_exact_impl.hpp"

#endif // guard
