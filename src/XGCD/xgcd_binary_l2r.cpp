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
void XGCD_BINARY_L2R(int64_t & G, int64_t & X, int64_t & Y, const int64_t & A, const int64_t & B){
  assert(X);
  assert(Y);

  const int64_t am = A >> 63;
  const int64_t bm = B >> 63;

  int64_t u1 = 1;
  int64_t u2 = 0;
  int64_t u3 = ANTL::negate_using_mask<int64_t>(am, A);
  int64_t v1 = 0;
  int64_t v2 = 1;
  int64_t v3 = ANTL::negate_using_mask<int64_t>(bm, B);
  
  // Swap u with v if u3 < v3.
  ANTL::cond_swap3_s64(u1, u2, u3, v1, v2, v3);
  while (v3 != 0) {
    int k = ANTL::msb_u64(u3) - ANTL::msb_u64(v3);

    // Subtract 2^k times v from u, and make sure u3 >= 0.
    uint64_t m;
   ANTL::sub_with_mask(m, u3, v3 << k);

    u1 -= v1 << k;
    u2 -= v2 << k;
    u1 = ANTL::negate_using_mask<int64_t>(m, u1);
    u2 = ANTL::negate_using_mask<int64_t>(m, u2);
    u3 = ANTL::negate_using_mask<int64_t>(m, u3);
    
    // Swap u with v if u3 < v3.
    ANTL::cond_swap3_s64(u1, u2, u3, v1, v2, v3);
  }

  if (u3 == ANTL::negate_using_mask<int64_t>(am,A)) {
    // a divides b.
    X = am | 1;  // either 1 or -1
    Y = 0;
  } else if (u3 == ANTL::negate_using_mask<int64_t>(bm, B)) {
    // b divides a.
    X = 0;
    Y = bm | 1;  // either 1 or -1
  } else {
#if (REDUCE_OUTPUT == 1)
    // Reduce u1 (mod b) and u2 (mod a) and correct for sign.
    int64_t q = u1 / b;
    X = ANTL::negate_using_mask<int64_t>(am, u1 - q*b);
    Y = ANTL::negate_using_mask<int64_t>(bm, u2 + q*a);
#else
    X = ANTL::negate_using_mask<int64_t>(am, u1);
    Y = ANTL::negate_using_mask<int64_t>(bm, u2);
#endif
  }
  G = u3;
}

template<>
void XGCD_BINARY_L2R_LEFT(int64_t & G, int64_t & X, const int64_t & A, const int64_t & B){
  assert(X);

  const int64_t am = A >> 63;
  const int64_t bm = B >> 63;

  int64_t u1 = 1;
  int64_t u3 = ANTL::negate_using_mask(am, A);
  int64_t v1 = 0;
  int64_t v3 = ANTL::negate_using_mask(bm, B);

  // Swap u with v if u3 < v3.
  ANTL::cond_swap2_s64(u1, u3, v1, v3);
  while (v3 != 0) {
    int k = ANTL::msb_u64(u3) - ANTL::msb_u64(v3);

    // Subtract 2^k times v from u, and make sure u3 >= 0.
    uint64_t m;
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
  assert(R2);
  assert(R1);
  assert(C2);
  assert(C1);
  assert(bound >= 0);
  int64_t r2 = R2;
  int64_t r1 = R1; 
  int64_t c2 = 0;
  int64_t c1 = -1;
  uint64_t s2 = r2 >> 63;
  uint64_t s1 = r1 >> 63;
  uint64_t cm = s2 ^ s1;
  r2 = ANTL::negate_using_mask(s2, r2);
  r1 = ANTL::negate_using_mask(s1, r1);
  assert(r2 >= r1);

  // Swap u with v if u3 < v3.
  Z ^= ANTL::cond_swap3_s64(c2, (int64_t &)s2, r2, c1, (int64_t &)s1, r1);
  while (r1 != 0 && r1 > bound) {
    int k = ANTL::msb_u64(r2) - ANTL::msb_u64(r1);

    // Subtract 2^k times r1 from r2, make sure r2 >= r1 >= 0
    uint64_t m;
    r2 = ANTL::sub_with_mask(m, r2, r1 << k);
    c2 -= ANTL::negate_using_mask(cm, c1 << k);

    r2 = ANTL::negate_using_mask(m, r2);
    s2 ^= m;
    cm ^= m;


    Z ^= ANTL::cond_swap3_s64(c2, (int64_t&)s2, r2, c1, (int64_t&)s1, r1);
  }

  R2 = ANTL::negate_using_mask(s2, r2);
  R1 = ANTL::negate_using_mask(s1, r1);
  C2 = c2;
  C1 = c1;
}

