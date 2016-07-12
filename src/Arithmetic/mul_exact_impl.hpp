/**
 * @file mul_exact_impl.hpp
 * @author Michael Jacobson
 * @remark template implementations of multiply functions computing only high-degree
 * terms.
 */

//
// PlainMulExact
//

template < class T >
void
PlainMulExact(T & x, const T & a, const T & b, long n)
{
  mul(x,a,b);
}



//
// PlainSqrExact
//

template < class T >
void
PlainSqrExact(T & x, const T & a, long n)
{
  sqr(x,a);
}



//
// MulExact
//

template < class T >
void
MulExact(T & x, const T & a, const T & b, long n)
{
  long size = deg(a) + deg(b);

#ifdef TRACE_MULEXACT
  cout << "MUL EXACT:  size=" << size << ", thresh=" << Thresholds<T>::get_mulexact_crossover(size,n) << ", n=" << n << endl;
#endif

  if (n == 0 || size < Thresholds<T>::get_mulexact_crossover(size,n)) {
    mul(x,a,b);
#ifdef TRACE_MULEXACT
    cout << "--> NTL" << endl;
#endif
  }
  else
    ::PlainMulExact(x,a,b,n);
}



//
// SqrExact
//

template < class T >
void
SqrExact(T & x, const T & a, long n)
{
  long size = deg(a) << 1;

#ifdef TRACE_MULEXACT
  cout << "SQR EXACT:  size=" << size << ", thresh=" << Thresholds<T>::get_sqrexact_crossover(size,n) << endl;
#endif

  if (n == 0 || size < Thresholds<T>::get_sqrexact_crossover(size,n)) {
    sqr(x,a);
#ifdef TRACE_MULEXACT
    cout << "--> NTL" << endl;
#endif
  }
  else
    ::PlainSqrExact(x,a,n);
}
