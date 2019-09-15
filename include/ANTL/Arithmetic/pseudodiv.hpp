/**
 * @file pseudodiv.hpp
 * @author Michael Jacobson
 * @brief Pseudodivision with polynomial base types
 */

#ifndef PSEUDODIV_H
#define PSEUDODIV_H

//#include <ANTL/utilities.hpp>
#include "../utilities.hpp"

// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.
NTL_CLIENT



//
// Pseudo division with remainder (for 1 inversion XGCD routines,
// polynomial base types only))
//

/**
 * @fn PseudoDivRem
 * Computes d in F, q, r in F[x] such that d a = q b + r and
 * deg(r) < deg(b).  Thus, q d^-1 and r d^-1 yield the usual
 * quotient and remainder.
 *
 * @param d The field element multiplier in the output
 * @param q The pseudo-quotient
 * @param r The pseudo-remainder
 * @param a The dividend (input)
 * @param b The divisor (input)
 */
template < class S, class T >
void PseudoDivRem(S & d, T & q, T & r, const T & a, const T & b);



/**
 * @fn PseudoDivRem_reduce
 * Computes d in F, q, r in F[x] such that d a = q b + r and
 * deg(r) < deg(b).  Thus, q d^-1 and r d^-1 yield the usual
 * quotient and remainder.
 *
 * In this version, coefficients that are multiplied by the leading
 * coefficient of b are reduced immediately.
 *
 * @param d The field element multiplier in the output
 * @param q The pseudo-quotient
 * @param r The pseudo-remainder
 * @param a The dividend (input)
 * @param b The divisor (input)
 */
template < class S, class T >
void PseudoDivRem_reduce(S & d, T & q, T & r, const T & a, const T & b);



/**
 * @fn PseudoDivRem_accumulate
 * Computes d in F, q, r in F[x] such that d a = q b + r and
 * deg(r) < deg(b).  Thus, q d^-1 and r d^-1 yield the usual
 * quotient and remainder.
 *
 * Unlike PseudoDivRem_reduce, this version does not immediately reduce coefficients
 * that are multiplied by the leading coefficient of b.  Instead, these are
 * allowed to accumulate and reduced once at the end of the function.
 *
 * @param d The field element multiplier in the output
 * @param q The pseudo-quotient
 * @param r The pseudo-remainder
 * @param a The dividend (input)
 * @param b The divisor (input)
 */
template < class S, class T >
void PseudoDivRem_accumulate(S & d, T & q, T & r, const T & a, const T & b);



//
// Declare specialized methods
//

template <> void PseudoDivRem(zz_p & d, zz_pX & q, zz_pX & r, const zz_pX & a, const zz_pX & b);
template <> void PseudoDivRem(ZZ_p & d, ZZ_pX & q, ZZ_pX & r, const ZZ_pX & a, const ZZ_pX & b);
template <> void PseudoDivRem(GF2E & d, GF2EX & q, GF2EX & r, const GF2EX & a, const GF2EX & b);

template <> void PseudoDivRem_reduce(zz_p & d, zz_pX & q, zz_pX & r, const zz_pX & a, const zz_pX & b);

template <> void PseudoDivRem_accumulate(zz_p & d, zz_pX & q, zz_pX & r, const zz_pX & a, const zz_pX & b);
template <> void PseudoDivRem_accumulate(ZZ_p & d, ZZ_pX & q, ZZ_pX & r, const ZZ_pX & a, const ZZ_pX & b);
template <> void PseudoDivRem_accumulate(GF2E & d, GF2EX & q, GF2EX & r, const GF2EX & a, const GF2EX & b);


// Unspecialized template definitions.
#include "../../../src/Arithmetic/pseudodiv_impl.hpp"

#endif // guard
