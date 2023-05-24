/**
 * @file xgcd_bin_l2r.hpp
 * @author Zack Baker
 * @brief XGCD methods using the Left-to-Right binary algorithms
 * 
 * @date 2022-05-04
 * 
 * 
 * 
 */

#ifndef XGCD_BIN_L2R_H
#define XGCD_BIN_L2R_H
#endif

#include <ANTL/common.hpp>

NTL_CLIENT

template <class T>
void XGCD_BINARY_L2R(T & G, T & X, T & Y, const T & A, const T & B);
template<>
void XGCD_BINARY_L2R(long & G, long & X, long & Y, const long & A, const long & B);

template <class T>
void XGCD_BINARY_L2R_LEFT(T & G, T & X, const T & A, const T & B);
template<>
void XGCD_BINARY_L2R_LEFT(long & G, long & X, const long & A, const long & B);

template <class T>
void XGCD_PARTIAL_BINARY_L2R(T & Z, T & R2, T & R1, T & C2, T & C1, const T bound);
template<>
void XGCD_PARTIAL_BINARY_L2R(long & Z, long & R2, long & R1, long & C2, long & C1, const long bound);
