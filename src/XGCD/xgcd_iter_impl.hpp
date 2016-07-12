/**
 * @file xgcd_iter_impl.hpp
 * @author Michael Jacobson
 * @remark template implementations of XGCD, XGCD_LEFT_ITER, XGCD_PARTIAL_ITER, and 
 * XGCD_PARTIAL_REDUCE_ITER
 */

//
// XGCD_ITER
//

template < class T >
void XGCD_ITER(T & G, T & X, T & Y, const T & A, const T & B)
{
  if (deg(B) < Thresholds<T>::get_pseudo_xgcd_crossover())
    XGCD_PSEUDO(G,X,Y,A,B);
  else
    XGCD_PLAIN(G,X,Y,A,B);
}



//
// XGCD_LEFT_ITER
//

template < class T >
void XGCD_LEFT_ITER(T & G, T & X, const T & A, const T & B)
{
  if (deg(B) < Thresholds<T>::get_pseudo_xgcd_left_crossover())
    XGCD_LEFT_PSEUDO(G,X,A,B);
  else
    XGCD_LEFT_PLAIN(G,X,A,B);
}



//
// Partial Euclidean algorithm (for NUCOMP)
//

template < class T >
void XGCD_PARTIAL_ITER(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag)
{  
  if (deg(R1) < Thresholds<T>::get_pseudo_xgcd_partial_crossover(deg(R1),bound))
    XGCD_PARTIAL_PSEUDO(R2,R1,C2,C1,bound,flag);
  else
    XGCD_PARTIAL_PLAIN(R2,R1,C2,C1,bound,flag);
}



//
// Partial Euclidean algorithm (for fast reduce) with revised termination for 
// polynomial types
//

template < class T >
void XGCD_PARTIAL_REDUCE_ITER(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even)
{  
  if (deg(R1) < Thresholds<T>::get_pseudo_xgcd_partial_crossover(deg(R1),bound))
    XGCD_PARTIAL_REDUCE_PSEUDO(R2,R1,B2,B1,bound,flag,even);
  else
    XGCD_PARTIAL_REDUCE_PLAIN(R2,R1,B2,B1,bound,flag,even);
}
