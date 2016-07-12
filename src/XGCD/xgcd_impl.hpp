/**
 * @file xgcd_impl.hpp
 * @author Michael Jacobson
 * @remark template implementations of XGCD, XGCD_LEFT, XGCD_PARTIAL, and 
 * XGCD_PARTIAL_REDUCE
 */

//
// XGCD_LEFT
//

template < class T >
void XGCD(T & G, T & X, T & Y, const T & A, const T & B)
{
#ifdef TRACE_XGCD
  cout << "XGCD:" << endl;
#endif
  if (deg(B) < Thresholds<T>::get_half_xgcd_crossover())
    XGCD_ITER(G,X,Y,A,B);
  else
    HXGCD(G,X,Y,A,B);
}



//
// XGCD_LEFT
//

template < class T >
void XGCD_LEFT(T & G, T & X, const T & A, const T & B)
{
#ifdef TRACE_XGCD
  cout << "XGCD_LEFT:" << endl;
#endif  
  if (deg(B) < Thresholds<T>::get_half_xgcd_left_crossover())
    XGCD_LEFT_ITER(G,X,A,B);
  else {
    cout << "Oh oh - trying to call HXGCD!" << endl;
    T Y;
    HXGCD(G,X,Y,A,B);
  }
}



//
// Partial Euclidean algorithm (for NUCOMP)
//

template < class T >
void XGCD_PARTIAL(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag)
{  
  /*
  if (deg(B) < Thresholds<T>::get_half_xgcd_partial_crossover())
    XGCD_PARTIAL_ITER(R2,R1,C2,C1,bound,flag);
  else
    HXGCD_PARTIAL(R2,R1,C2,C1,bound,flag);
  */

  XGCD_PARTIAL_ITER(R2,R1,C2,C1,bound,flag);
}



//
// Partial Euclidean algorithm (for fast reduce) with revised termination for 
// polynomial types
//

template < class T >
void XGCD_PARTIAL_REDUCE(T & R2, T & R1, T & B2, T & B1, long bound, bool & flag, bool even)
{  
  /*
  if (deg(B) > Thresholds<T>::get_half_xgcd_partial_crossover())
    XGCD_PARTIAL_REDUCE_ITER(R2,R1,B2,B1,bound,flag,even);
  else
    HXGCD_PARTIAL_REDUCE(R2,R1,B2,B1,bound,flag,even);
  */

  XGCD_PARTIAL_REDUCE_ITER(R2,R1,B2,B1,bound,flag,even);
}
