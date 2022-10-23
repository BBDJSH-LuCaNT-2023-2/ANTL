/**
 * @file common.hpp
 * @author Michael Jacobson
 * @brief General-purpose methods to interface with NTL.  All library files should include this file.
 */

#ifndef ANTL_COMMON_H
#define ANTL_COMMON_H

#include <cmath>
#include <gmp.h>

#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/GF2EX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZ_limbs.h>

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>

#include <boost/math/bindings/rr.hpp>
using boost::math::ntl::atan;

// We use the NTL namespace everywhere. Rather than have a using directive in
// every file, we just put it here, for convenience and clarity.
NTL_CLIENT


/**
 * @brief Suggested beginning of the definition of a function-like macro.
 *
 * This helps to prevent heisenbugs caused by macros that have tricky semantics
 * when placed in, say, an if-block. The optimizer strips out the do-while.
 */
#define ANTL_BEGIN_MACRO_FUNCTION	do {

/**
 * @brief Used to close a macro definition started with \c ANTL_BEGIN_MACRO.
 */
#define ANTL_END_MACRO_FUNCTION		} while(0)

// Include macros for debugging output, and general output.
#include <ANTL/debug.hpp>

/**
 * @brief Null macro used to mark a by-reference argument pass.
 *
 * This macro improves human readability by indicating that a method call
 * may change an argument in the caller scope.
 */
#define REF

// Forward declarations for prototypes in this file.
namespace NTL {
  //
  // NTL-like methods, for completeness. There exists a "set" and "clear" for most (all?) NTL types.
  // The below definitions extend this functionality for other types.

  inline void clear (long &X)        { X = 0; }
  inline void set (long &X)          { X = 1; }
  inline void clear (float &i)       { i = 0.0f; }
  inline void set (float &i)         { i = 1.0f; }
  inline void clear (double &i)      { i = double(0); }
  inline void set (double &i)        { i = double(1); }
  inline void clear (quad_float & i) { i = double(0); }
  inline void set (quad_float & i)   { i = double(1); }


  inline long IsOne (const long &X)         { return (X == 1); }
  inline long IsZero (const long &X)        { return (X == 0); }
  inline long IsOdd (const long &X)         { return (X & 1); }
  inline long sign (const long &X)          { return (X > 0) ? 1 : ((X < 0) ? -1 : 0); }

  inline long IsOne (const float &X)         { return (X == 1.0f); }
  inline long IsZero (const float &X)        { return (X == 0.0f); }

  inline long IsOne (const double &X)         { return (X == double(1)); }
  inline long IsZero (const double &X)        { return (X == double(0)); }

  inline long IsOne (const quad_float &X)         { return (X == double(1)); }
  inline long IsZero (const quad_float &X)        { return (X == double(0)); }

  // JLM: needed to compile quadratic_form. I have not tested them
  // for correctness.
  inline long deg(long X)              { return 0; }
  inline long deg(const ZZ & X)        { return 0; }
  inline long LeadCoeff(long X)        { return X; }
  inline ZZ LeadCoeff(const ZZ & X)    { return X; }
  inline void MakeMonic(long & X)      { X = 1; }
  inline void MakeMonic(ZZ & X)        { X = to_ZZ(1); }

  // assign(C,A) method for NTL classes that only support the assignment operator
  // (for compatibility with templated classes like exponentiation)
  template < class T >
  inline void assign(T &C, const T &A) { C = A; }

  inline void SqrRoot(double &x, const double &a) { x = std::sqrt(a); }
  inline void ComputePi(double & p){ p = 3.141592653589793;}
  // procedural arithmetic operations for standard types, to increase compatibility with NTL
  inline void add ( long &C, const long &A, const long &B )     { C = A + B; }
  inline void add ( float &C, const float &A, const float &B )     { C = A + B; }
  inline void add ( double &C, const double &A, const double &B )     { C = A + B; }

  inline void sub ( long &C, const long &A, const long &B )     { C = A - B; }
  inline void sub ( float &C, const float &A, const float &B )     { C = A - B; }
  inline void sub ( double &C, const double &A, const double &B )     { C = A - B; }

  inline void mul ( long &C, const long &A, const long &B )     { C = A * B; }
  inline void mul ( float &C, const float &A, const float &B )     { C = A * B; }
  inline void mul ( double &C, const double &A, const double &B )     { C = A * B; }

  inline void div ( long &C, const long &A, const long &B )     { C = A / B; }
  inline void div ( float &C, const float &A, const float &B )     { C = A / B; }
  inline void div ( double &C, const double &A, const double &B )     { C = A / B; }

  inline void sqr ( long &C, const long &A)     { C = A * A; }
  inline void sqr ( float &C, const float &A)     { C = A * A; }
  inline void sqr ( double &C, const double &A)     { C = A * A; }

  inline void power(double &C, double &A, double &E) { C= std::pow(A,E); }
  inline void power(float &C, float &A, float &E) { C= std::pow(A,E); }

  inline void pow(double &C, double &A, const double &E) { C= std::pow(A,E); }
  inline void pow(float &C, float &A, const float &E) { C= std::pow(A,E); }


  inline void abs(long & z, const long a) { z =  std::abs(a); }
  inline void abs(double & z, const double a) { z = std::abs(a); }

  inline void floor(long & z, const long a) { z =  std::floor(a); }
  inline void floor(double & z, const double a) { z = std::floor(a); }

  inline void log(double & z, const double a) { z = std::log(a); }

  inline void atan_val(RR & z, const RR & a){z = atan(a).value(); }
  inline void atan_val(double & z, const double & a){z = atan(a); }
}

namespace ANTL {
  inline int SqrRoot (const int & a) { return (int) ::floor(::sqrt((double) a)); }

  void DivRem (long &q, long &r, long a, long b);
  inline long SqrRoot (const long &a) { return (long) ::floor(::sqrt((double) a)); }

  //
  // cmath-like methods, for completeness.
  // Most of these exist because of problems where gcc doesn't find the
  // cmath methods, presumably because of namespace issues).
  //
  inline float sqrt(float x) { return (float)std::sqrt((double)x); }
  inline double sqrt(double x) { return std::sqrt(x); }
  inline long abs(long x) { return std::abs(x); }
  inline float abs(float x) { return (float)std::abs((double)x); }
  inline double abs(double x) { return std::abs(x); }
  inline quad_float abs(const quad_float & x) { return NTL::fabs(x); }
  inline float exp(float x) { return (float)std::exp((double) x); }
  inline double exp(double x) { return std::exp(x); }
  inline float log(float x) { return (float)std::log((double)x); }
  inline double log(double x) { return std::log(x); }
  inline float pow(float x, float e) { return (float)std::pow((double)x, (double)e ); }
  inline double pow(double x, float e) { return std::pow(x,e); }

  //
  // quadratic residuosity functions
  //

  /* Jacobi functions - assumes a is reduced mod n */
  long Jacobi_base (const long & a, const long & n);
  long Jacobi_base (const long long & a, const long long & n);
  long Jacobi_base (const ZZ & a, const ZZ & n);

  /* Jacobi functions - no preconditions on a and n */
  long Jacobi(const long & a, const long & n);
  long Jacobi(const long long & a, const long long & n);
  long Jacobi(const long long & a, const long & n);
  long Jacobi(const ZZ & a, const ZZ & n);
  long Jacobi(const ZZ & a, const long & n);
  long Jacobi(const ZZ_pX & a, const ZZ_pX & n);
  long Jacobi(const zz_pX & a, const zz_pX & n);
  long Jacobi(const ZZ_pEX & a, const ZZ_pEX & n);
  long Jacobi(const zz_pEX & a, const zz_pEX & n);
//   long Jacobi(const GF2EX & h, const GF2EX & a, const GF2EX & n);

  /* Kronecker functions - no preconditions on a and n */
  long Kronecker(const long & a, const long & n);
  long Kronecker(const long long & a, const long long & n);
  long Kronecker(const long long & a, const long & n);
  long Kronecker(const ZZ & a, const ZZ & n);
  long Kronecker(const ZZ & a, const long & n);
  long Kronecker(const ZZ_pX & a, const ZZ_pX & n);
  long Kronecker(const zz_pX & a, const zz_pX & n);
  long Kronecker(const ZZ_pEX & a, const ZZ_pEX & n);
  long Kronecker(const zz_pEX & a, const zz_pEX & n);
  long Kronecker(const GF2EX & h, const GF2EX & f, const GF2EX & n);

  //mask negation functions: m must be 0 or -1
  template <class T>
  T negate_using_mask(const uint64_t m, const T x);

  int64_t sub_with_mask(uint64_t & m, const int64_t & a, const int64_t & b);

  void cond_swap2_s64(int64_t & u1, int64_t & u2, int64_t & v1, int64_t & v2);
  uint64_t cond_swap3_s64(int64_t & u1,
				      int64_t & u2,
				      int64_t & u3,
				      int64_t & v1,
				      int64_t & v2,
				      int64_t & v3);

  int msb_u64(uint64_t x);
  // finite field cardinality macros
  template <class> ZZ CARDINALITY(void);

  template <> inline ZZ CARDINALITY < ZZ > (void) {
    return to_ZZ(0);
  }

  template <> inline ZZ CARDINALITY < long > (void) {
    return to_ZZ(0);
  }

  template <> inline ZZ CARDINALITY < ZZ_p > (void) {
    return ZZ_p::modulus ();
  }

  template <> inline ZZ CARDINALITY < zz_p > (void) {
    return to_ZZ (zz_p::modulus ());
  }

  template <> inline ZZ CARDINALITY < ZZ_pX > (void) {
    return ZZ_p::modulus ();
  }

  template <> inline ZZ CARDINALITY < zz_pX > (void) {
    return to_ZZ (zz_p::modulus ());
  }

  template <> inline ZZ CARDINALITY < ZZ_pE > (void) {
    return ZZ_pE::cardinality ();
  }

  template <> inline ZZ CARDINALITY < zz_pE > (void) {
    return zz_pE::cardinality ();
  }

  template <> inline ZZ CARDINALITY < ZZ_pEX > (void) {
    return ZZ_pE::cardinality ();
  }

  template <> inline ZZ CARDINALITY < zz_pEX > (void) {
    return zz_pE::cardinality ();
  }

  template <> inline ZZ CARDINALITY < GF2E > (void) {
    return GF2E::cardinality ();
  }

  template <> inline ZZ CARDINALITY < GF2EX > (void) {
    return GF2E::cardinality ();
  }



  //
  // Template-friendly methods for common conversions & constants.
  // e.g.: ZZ x = to<ZZ>(myLongVariable);
  //


  template <class T> T to(const int& a)          { return T(a); }
  template <class T> T to(const long& a)         { return T(a); }
  template <class T> T to(const float& a)        { return T(a); }
  template <class T> T to(const double& a)       { return T(a); }
  template <class T> T to(const ZZ& a)           { return T(a); }
  template <class T> T to(const RR& a)           { return T(a); }
  template <class T> T to(const quad_float & a)  { return T(a); }
  template <class T> T to(const GF2EX & a)       { return T(a); }
  template <class T> T to(const zz_pX & a)       { return T(a); }
  template <class T> T to(const ZZ_pX & a)       { return T(a); }
  template <class T> T to(const zz_pEX & a)      { return T(a); }
  template <class T> T to(const ZZ_pEX & a)      { return T(a); }

  template<> inline int to<int>(const int& a)        { return a; }
  template<> inline int to<int>(const long& a)       { return to_int(a); }
  template<> inline int to<int>(const float& a)      { return to_int(a); }
  template<> inline int to<int>(const double& a)     { return to_int(a); }
  template<> inline int to<int>(const ZZ& a)         { return to_int(a); }
  template<> inline int to<int>(const RR& a)         { return to_int(a); }
  template<> inline int to<int>(const quad_float& a) { return to_int(a); }

  template<> inline long to<long>(const int& a)        { return to_long(a); }
  template<> inline long to<long>(const long& a)       { return a; }
  template<> inline long to<long>(const float& a)      { return to_long(a); }
  template<> inline long to<long>(const double& a)     { return to_long(a); }
  template<> inline long to<long>(const ZZ& a)         { return to_long(a); }
  template<> inline long to<long>(const RR& a)         { return to_long(a); }
  template<> inline long to<long>(const quad_float& a) { return to_long(a); }

  template<> inline float to<float>(const int& a)        { return to_float(a); }
  template<> inline float to<float>(const long& a)       { return to_float(a); }
  template<> inline float to<float>(const float& a)      { return a; }
  template<> inline float to<float>(const double& a)     { return to_float(a); }
  template<> inline float to<float>(const ZZ& a)         { return to_float(a); }
  template<> inline float to<float>(const RR& a)         { return to_float(a); }
  template<> inline float to<float>(const quad_float& a) { return to_float(a); }

  template<> inline double to<double>(const int& a)        { return to_double(a); }
  template<> inline double to<double>(const long& a)       { return to_double(a); }
  template<> inline double to<double>(const float& a)      { return to_double(a); }
  template<> inline double to<double>(const double& a)     { return a; }
  template<> inline double to<double>(const ZZ& a)         { return to_double(a); }
  template<> inline double to<double>(const RR& a)         { return to_double(a); }
  template<> inline double to<double>(const quad_float& a) { return to_double(a); }

  template<> inline ZZ to<ZZ>(const int& a)        { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const long& a)       { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const float& a)      { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const double& a)     { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const ZZ& a)         { return a; }
  template<> inline ZZ to<ZZ>(const RR& a)         { return to_ZZ(a); }
  template<> inline ZZ to<ZZ>(const quad_float& a) { return to_ZZ(a); }

  template<> inline RR to<RR>(const int& a)        { return to_RR(a); }
  template<> inline RR to<RR>(const long& a)       { return to_RR(a); }
  template<> inline RR to<RR>(const float& a)      { return to_RR(a); }
  template<> inline RR to<RR>(const double& a)     { return to_RR(a); }
  template<> inline RR to<RR>(const ZZ& a)         { return to_RR(a); }
  template<> inline RR to<RR>(const RR& a)         { return a; }
  template<> inline RR to<RR>(const quad_float& a) { return to_RR(a); }

  template<> inline quad_float to<quad_float>(const int& a)        { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const long& a)       { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const float& a)      { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const double& a)     { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const ZZ& a)         { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const RR& a)         { return to_quad_float(a); }
  template<> inline quad_float to<quad_float>(const quad_float& a) { return a; }

  template<> inline GF2EX to<GF2EX>(const int& a)        { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const long& a)       { return to_GF2EX(a); }
  template<> inline GF2EX to<GF2EX>(const float& a)      { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const double& a)     { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const ZZ& a)         { return to_GF2EX(a); }
  template<> inline GF2EX to<GF2EX>(const RR& a)         { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const quad_float& a) { return to_GF2EX(to_ZZ(a)); }
  template<> inline GF2EX to<GF2EX>(const GF2EX& a) { return a; }

  template<> inline ZZ_pX to<ZZ_pX>(const ZZ_pX& a) { return a; }
  template<> inline ZZ_pEX to<ZZ_pEX>(const ZZ_pEX& a) { return a; }
  template<> inline zz_pX to<zz_pX>(const zz_pX& a) { return a; }
  template<> inline zz_pEX to<zz_pEX>(const zz_pEX& a) { return a; }

  template<> inline zz_pX to<zz_pX>(const int & a) { return zz_pX(a,0); }

  void ZZToMpz(const ZZ & A, mpz_t & a);
  void MpzToZZ(const mpz_t & a, ZZ & A);
} // ANTL



//
// Utility routines for polynomials.
//

/**
 * @fn get_poly_modq
 * Computes a polynomial mod q corresponding to the integer X (q-adic
 * representation of X).
 *
 * @param p Reference to an instance of GF2EX to hold the newly
 * created polynomial.
 * @param X The integer to create the polynomial from.
 * @param q The modulus of the polynomial creation.
 */
void get_poly_modq (GF2X & p, const ZZ & X, const ZZ & q);
void get_poly_modq (GF2EX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (ZZ_pX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (ZZ_pEX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (zz_pX & p, const ZZ & X, const ZZ & q);
void get_poly_modq (zz_pEX & p, const ZZ & X, const ZZ & q);


/**
 * @brief Computes the integer X corresponding to the q-adic
 * representation of p.
 */
template < class T >
ZZ eval_poly (const T & p, const ZZ & q);
template <> ZZ eval_poly <ZZ_pEX> (const ZZ_pEX & p, const ZZ & q);
template <> ZZ eval_poly <zz_pEX> (const zz_pEX & p, const ZZ & q);
template <> ZZ eval_poly <GF2EX> (const GF2EX & p, const ZZ & q);



//
// quadratic residuosity functions
//

/* Jacobi functions - assumes a is reduced mod n */
long Jacobi_base (const ZZ & a, const ZZ & n);
long Jacobi_base (const long & a, const long & n);

/* Jacobi functions - no preconditions on a and n */
//long Jacobi(const ZZ & a, const ZZ & n);
//long Jacobi(const long & a, const long & n) { return long(Jacobi( to_ZZ(a), to_ZZ(n))); }
long Jacobi(const ZZ_pX & a, const ZZ_pX & n);
long Jacobi(const zz_pX & a, const zz_pX & n);
long Jacobi(const ZZ_pEX & a, const ZZ_pEX & n);
long Jacobi(const zz_pEX & a, const zz_pEX & n);

/* GF2X quadratic character */
/**
 * @brief Computes the trace of a (mod p), i.e.,
 * 		tr = sum_{i=0}^{vt-1} a^2^i (mod p)
 *  where deg(p) = v.
 */
void trace (GF2X & t, const GF2X & a, const GF2X & pmod);
long Jacobi (const GF2X & h, const GF2X & f, const GF2X & n);

/* GF2EX quadratic character */
/**
 * @brief Computes the trace of a (mod p), i.e.,
 * 		tr = sum_{i=0}^{vt-1} a^2^i (mod p)
 *  where deg(p) = v.
 */
void trace (GF2EX & t, const GF2EX & a, const GF2EXModulus & pmod);
long Jacobi (const GF2EX & h, const GF2EX & f, const GF2EX & n);



//
// modular square root functions
//

//
// get_qdp (computes a square root of Delta mod p)
//   - used by templated ressol)
//
void get_qdp(ZZ & q, const ZZ_pX & p);
void get_qdp(ZZ & q, const zz_pX & p);
void get_qdp(ZZ & q, const ZZ_pEX & p);
void get_qdp(ZZ & q, const zz_pEX & p);

/**
 * @brief Computes x, a a solution of x^2 = a (mod p).
 * @return
 * 		 0  if the equation has a unique solution (p | a)
 * 		-1  if x^2 = a (mod p) has no solutions
 * 		 1  if x^2 = a (mod p) has 2 solutions
 */
template < class T >
long ressol (T & x, const T & a, const T & p);



/**
 * @brief Computes x, a a solution of x^2 + hx -f = 0 (mod p).
 * @return
 * 		 0  if the equation has a unique solution (p | h)
 * 		-1  if x^2 + hx -f = 0 (mod p) has no solutions
 * 		 1  if x^2 + hx -f = 0 (mod p) has 2 solutions
 */
long ressol (GF2EX & x, const GF2EX & h, const GF2EX & f, const GF2EX & p);


/**
* @brief Checks whether all coordinates of vec1 and vec2 are within maxdist of each other
* Currently for RR and doubles
* @param[out] true if abs(v1[i] -v2[i]) < maxdist for all i
* @param[in] B vec1 is a vector of real type entries.
* @param[in] vec2 is a vector with real type entries
* @param[in] maxdist is a real type value that indicates the max acceptable distance of entries to be considered close
* @pre vec1 and vec2 need to be the same size.
*/
template<typename PType>
bool is_close(const std::vector<PType> & vec1, const std::vector<PType> & vec2, const PType & max_dist);


/*
* @brief used in the computation of compact representations. Not yet implemented
*/
template<typename Type, typename PType>
void compute_initial_s(const std::vector<PType> & alpha, const int k_bound);

/*
* @brief Takes a length r log vector and returns the corresponding length r+1 exponentiated vector. equivalent to create_target from the pari stuff.
* @param[in] log_vec is a vector with real entries having length equal to the unit rank r of the field
* @param[out] valuationvec is a length r+1 vector, consisting of  exp(log_vec[i]), -sum(deg(i)*log_vec[i])/deg(r+1) )
* Here deg(i) is 1 if the ith coordinate is a real embedding, 2 otherwise.
*/
template<typename PType>
void log_to_valuation(std::vector<PType> &valuation_vec, const std::vector<PType> & log_vec, const int r1);



// Unspecialized template definitions.
#include "../../src/common_impl.hpp"


#endif // guard
