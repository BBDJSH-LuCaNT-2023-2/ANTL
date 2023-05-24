/**
 * @file xgcd_binary_l2r.cpp
 * @author Zack Baker
 * @brief Specializations of left-to-right binary xgcd algorithms
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <ANTL/XGCD/xgcd_binary_l2r.hpp>
#include <ANTL/common.hpp>


template<>
void XGCD_BINARY_L2R(long & G, long & X, long & Y, const long & A, const long & B){
//   std::cout << "using XGCD_BINARY_L2R" << std::endl;


  const long am = A >> 63;
  const long bm = B >> 63;

  long u1 = 1;
  long u2 = 0;
  long u3 = ANTL::negate_using_mask(am, A);
  long v1 = 0;
  long v2 = 1;
  long v3 = ANTL::negate_using_mask(bm, B);
  
  // Swap u with v if u3 < v3.
  ANTL::cond_swap3_s64(u1, u2, u3, v1, v2, v3);
  while (v3 != 0) {
    int k = ANTL::msb_u64(u3) - ANTL::msb_u64(v3);

    // Subtract 2^k times v from u, and make sure u3 >= 0.
    ulong m;
    u3 = ANTL::sub_with_mask(m, u3, v3 << k);
    u1 -= v1 << k;
    u2 -= v2 << k;
    u1 = ANTL::negate_using_mask(m, u1);
    u2 = ANTL::negate_using_mask(m, u2);
    u3 = ANTL::negate_using_mask(m, u3);
    
    // Swap u with v if u3 < v3.
    ANTL::cond_swap3_s64(u1, u2, u3, v1, v2, v3);
  }

  if (u3 == ANTL::negate_using_mask(am, A)) {
    // a divides b.
    X = am | 1;  // either 1 or -1
    Y = 0;
  } else if (u3 == ANTL::negate_using_mask(bm, B)) {
    // b divides a.
    X = 0;
    Y = bm | 1;  // either 1 or -1
  } else {
#if (REDUCE_OUTPUT == 1)
    // Reduce u1 (mod b) and u2 (mod a) and correct for sign.
    long q = u1 / b;
    X = ANTL::negate_using_mask_s64(am, u1 - q*b);
    Y = ANTL::negate_using_mask_s64(bm, u2 + q*a);
#else
    X = ANTL::negate_using_mask(am, u1);
    Y = ANTL::negate_using_mask(bm, u2);
#endif
  }
  
  G=u3;
}


template<>
void XGCD_BINARY_L2R_LEFT(long & G, long & X, const long & A, const long & B){
//   std::cout << "using XGCD_BINARY_L2R_LEFT" << std::endl;
//   assert(X);

  const long am = A >> 63;
  const long bm = B >> 63;

  long u1 = 1;
  long u3 = ANTL::negate_using_mask(am, A);
  long v1 = 0;
  long v3 = ANTL::negate_using_mask(bm, B);

  // Swap u with v if u3 < v3.
  ANTL::cond_swap2_s64(u1, u3, v1, v3);
  while (v3 != 0) {
    int k = ANTL::msb_u64(u3) - ANTL::msb_u64(v3);

    // Subtract 2^k times v from u, and make sure u3 >= 0.
    ulong m;
    u3 = ANTL::sub_with_mask(m, u3, v3 << k);
    u1 -= v1 << k;
    u1 = ANTL::negate_using_mask(m, u1);
    u3 = ANTL::negate_using_mask(m, u3);

    // Swap u with v if u3 < v3.
    ANTL::cond_swap2_s64(u1, u3, v1, v3);
  }

  if (u3 == ANTL::negate_using_mask(am, A)) {
    // a divides b.
    X = am | 1;  // either 1 or -1
  } else if (u3 == ANTL::negate_using_mask(bm, B)) {
    // b divides a.
    X = 0;
  } else {
    // Reduce u1 (mod b/u3) and correct for sign.
    X = ANTL::negate_using_mask(am, u1 % B);
  } 
  G = u3;
}

template<>
void XGCD_PARTIAL_BINARY_L2R(long & Z, long & R2, long & R1, long & C2, long & C1, const long bound){
//     std::cout << "using XGCD_PARTIAL_BINARY_L2R" << std::endl;

//   assert(bound >= 0);
  int k = 0;
  long R2_orig = R2;
  long r2 = R2;
  long r1 = R1;
  long c2 = 0;
  long c1 = -1;
  long s2 = r2 >> 63;
  long s1 = r1 >> 63;
  ulong cm = s2 ^ s1, m;
  Z = 0;

//   assert(r2 >= r1);

  // Swap u with v if u3 < v3.
  Z ^= ANTL::cond_swap3_s64(c2, s2, r2, c1, s1, r1);

  while (r1 > bound) {
    k = ANTL::msb_u64(r2) - ANTL::msb_u64(r1);

    // Subtract 2^k times r1 from r2, make sure r2 >= r1 >= 0
    m = 0;
    r2 = ANTL::sub_with_mask(m, r2, r1 << k);
    c2 -= c1 << k;

    while(r2 < 0) {
      r2 += m & r1;
      c2 += m & c1;
    }

    Z ^= ANTL::cond_swap3_s64(c2, (long&)s2, r2, c1, (long&)s1, r1);
  }
  R2 = r2;
  R1 = r1;
  C2 = c2;
  C1 = c1;

//   assert(abs(R1*C2 - R2*C1) == R2_orig);

}

