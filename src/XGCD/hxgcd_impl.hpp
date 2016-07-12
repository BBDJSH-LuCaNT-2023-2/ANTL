/**
 * @file hxgcd_impl.hpp
 * @author Laurent Imbert
 * @brief
 */



template <class T>
class TMatrix {
private:

   TMatrix(const TMatrix&);  // disable
   T elts[2][2];

public:

   TMatrix() { }
   ~TMatrix() { }

   void operator=(const TMatrix&);
   T& operator() (long i, long j) { return elts[i][j]; }
   const T& operator() (long i, long j) const { return elts[i][j]; }
};

template <class T>
void TMatrix<T>::operator=(const TMatrix& M)
{
   elts[0][0] = M.elts[0][0];
   elts[0][1] = M.elts[0][1];
   elts[1][0] = M.elts[1][0];
   elts[1][1] = M.elts[1][1];
}

// -----------------------------------------------------------------------------

template <class T>
void mul(T& A, T& B, const TMatrix<T>& M)
// (A, B) = M*(A, B)^T
{
   T R1, R2;

   mul(R1, M(0,0), A);
   mul(R2, M(0,1), B);
   add(R2, R1, R2);

   mul(R1, M(1,0), A);
   A = R2;
   mul(R2, M(1,1), B);
   add(B, R1, R2);

   R1.kill();
   R2.kill();
}

template <class T>
void mul(TMatrix<T>& A, TMatrix<T>& B, TMatrix<T>& C)
// A = B*C, B and C are destroyed
{
   T B00, B01, B10, B11, C0, C1, T1, T2;

   mul(T1, B(0,0), C(0,0));
   mul(T2, B(0,1), C(1,0));
   add(A(0,0), T1, T2);

   mul(T1, B(1,0), C(0,0));
   mul(T2, B(1,1), C(1,0));
   add(A(1,0), T1, T2);

   C(0,0).kill();
   C(1,0).kill();

   mul(T1, B(0,0), C(0,1));
   mul(T2, B(0,1), C(1,1));
   add(A(0,1), T1, T2);

   mul(T1, B(1,0), C(0,1));
   mul(T2, B(1,1), C(1,1));
   add(A(1,1), T1, T2);

   C(0,1).kill();
   C(1,1).kill();

   B(0,0).kill();
   B(1,0).kill();
   B(0,1).kill();
   B(1,1).kill();
}

template <class T>
void strassen_mul (TMatrix<T>& C, TMatrix<T>& A, TMatrix<T>& B)
// C = A*B, A and B are destroyed
/* we follow the code from SAGE 1.6, file strassen.pyx */
{
  T S0, T0, S1, T1, S2, T2, S3, T3, P0, P1, P2, P3, P4, P5, P6;

  add (S0, A(1,0), A(1,1));
  add (T0, B(0,1), B(0,0));
  add (S1, S0, A(0,0));
  add (T1, B(1,1), T0);
  add (S2, A(0,0), A(1,0));
  add (T2, B(1,1), B(0,1));
  add (S3, A(0,1), S1);
  add (T3, B(1,0), T1);
  mul (P0, A(0,0), B(0,0));
  mul (P1, A(0,1), B(1,0));
  mul (P2, S0, T0);
  mul (P3, S1, T1);
  mul (P4, S2, T2);
  mul (P5, S3, B(1,1));
  mul (P6, A(1,1), T3);
  add (C(0,0), P0, P1);     /* U0 */
  add (C(0,1), P0, P3);     /* U1 */
  add (C(1,1), C(0,1), P4); /* U2 */
  add (C(1,0), C(1,1), P6); /* U3 */
  add (C(1,1), C(1,1), P2); /* U4 */
  add (C(0,1), C(0,1), P2); /* U5 */
  add (C(0,1), C(0,1), P5); /* U6 */
  A(0,0).kill();
  A(1,0).kill();
  A(0,1).kill();
  A(1,1).kill();
  B(0,0).kill();
  B(1,0).kill();
  B(0,1).kill();
  B(1,1).kill();
}


// -----------------------------------------------------------------------------

// Computes iteratively a 2 x 2 polynomial matrix M_out such that
// M_out * (A,B)^T = (A',B')^T, where A', B' are consecutive polynomials
// in the Euclidean remainder sequence of A,B , and B' is the polynomial
// of highest degree satisfying deg(B') <= deg(A) - d_red.
// Replace (A,B) with (A',B')
template <class T>
void Iter_HALF_GCD(TMatrix<T>& M_out, T& A, T& B, long d_red)
{
  M_out(0,0).SetMaxLength(d_red);
  M_out(0,1).SetMaxLength(d_red);
  M_out(1,0).SetMaxLength(d_red);
  M_out(1,1).SetMaxLength(d_red);

  set(M_out(0,0));   clear(M_out(0,1));
  clear(M_out(1,0)); set(M_out(1,1));
  
  long goal = deg(A) - d_red;
  
  if (deg(B) <= goal)
    return;
  
  T Q, t(INIT_SIZE, d_red);

  while (deg(B) > goal)
    {
      DivRem(Q, A, A, B);
      swap(A, B);

      mul(t, Q, M_out(1,0)); // t = Q.v
      sub(t, M_out(0,0), t); // t = u - Q.v
      M_out(0,0) = M_out(1,0);  // u = v
      M_out(1,0) = t; // v = t

      mul(t, Q, M_out(1,1));
      sub(t, M_out(0,1), t);
      M_out(0,1) = M_out(1,1);
      M_out(1,1) = t;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// New version

template < class T >
void HALF_GCD(TMatrix<T>& M_out, T& A, T& B, long d_red)
{
   long degA = deg(A);   

   if (IsZero(B) || deg(B) <= degA - d_red)
     {
       set(M_out(0,0));   clear(M_out(0,1));
       clear(M_out(1,0)); set(M_out(1,1));
       return;
     }

   if (d_red <= Thresholds<T>::hxgcd_inner_crossover)
     {
       Iter_HALF_GCD(M_out, A, B, d_red);
       return;
     }

   long d1 = (d_red + 1) / 2;       /* d1 is about m/4 */
   if (d1 < 1) d1 = 1;
   if (d1 >= d_red) d1 = d_red - 1;

   long n = degA - 2 * d1 + 2;      /* n is the ignored part, about m/2 */
   if (n < 0) n = 0;

   T A1, B1;

   if (n != 0)
     {
       trunc(A1, A, n);          // A1 = A mod x^n
       trunc(B1, B, n);          // B1 = B mod x^n
       RightShift (A, A, n);     // A = A div x^n has about m/2 bits
       RightShift (B, B, n);     // B = B div x^n has about m/2 bits
     }

   TMatrix<T> M1;

   HALF_GCD (M1, A, B, d1);
   // the entries of M1 have m/4 bits, and A,B have been reduced to m/4 bits

   if (n != 0)
     {
       LeftShift (A, A, n);
       LeftShift (B, B, n);
       mul (A1, B1, M1);       // A1,B1: m/2  M1:m/4  cost: 4 M(m/2,m/4)
       add (A, A, A1);
       add (B, B, B1);
     }

   /* now A, B have 3m/4 bits */

   long d2 = deg(B) - (degA - d_red); /* should be about m/2 */

   if (IsZero(B) || d2 <= 0)
     {
       M_out = M1;
       return;
     }

   T Q;
   TMatrix<T> M2;

   DivRem(Q, A, A, B);
   swap(A, B);

   T t(INIT_SIZE, deg(M1(1,1))+deg(Q)+1);

   mul(t, Q, M1(1,0));
   sub(t, M1(0,0), t);
   swap(M1(0,0), M1(1,0));
   swap(M1(1,0), t);

   t.kill();

   t.SetMaxLength(deg(M1(1,1))+deg(Q)+1);

   mul(t, Q, M1(1,1));
   sub(t, M1(0,1), t);
   swap (M1(0,1), M1(1,1));
   swap (M1(1,1), t);

   if (IsZero(B) || deg(B) <= degA - d_red)
     {
       M_out = M1;
       return;
     }

   n = deg(A) - 2 * d2 + 2;     // should be about m/4
   if (n < 0) n = 0;

   if (n != 0)
     {
       trunc(A1, A, n);
       trunc(B1, B, n);        // A1,B1 have m/4 bits
       RightShift(A, A, n);    // A,B have m/2 bits
       RightShift(B, B, n);
     }

   HALF_GCD (M2, A, B, d2);
   // now A,B have m/4 bits, like the entries of M2

   if (n != 0)
     {
       LeftShift (A, A, n);
       LeftShift (B, B, n);
       mul (A1, B1, M2);
       add (A, A, A1);
       add (B, B, B1);
     }

   //   mul (M_out, M2, M1); 
   strassen_mul(M_out, M2, M1);
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


// Computes recursively a 2 x 2 polynomial matrix M_out such that
// M_out * (A,B)^T = (A',B')^T, where A', B' are consecutive polynomials
// in the Euclidean remainder sequence of A,B , and B' is the polynomial
// of highest degree satisfying deg(B') <= deg(A) - d_red.
// Replace (A,B) with (A',B')
// Usually if deg(A) = m, then d_red = m/2
// The crossover parameter is for profiling purpose only!


// template < class T >
// void HALF_GCD(TMatrix<T>& M_out, T& A, T& B, long d_red, long crossover=HalfGCD_Crossover<T>())
// {
//    long degA = deg(A);   

//    if (IsZero(B) || deg(B) <= degA - d_red)
//      {
//        set(M_out(0,0));   clear(M_out(0,1));
//        clear(M_out(1,0)); set(M_out(1,1));
//        return;
//      }

//    if (d_red <= crossover)
//      {
//        Iter_HALF_GCD(M_out, A, B, d_red);
//        return;
//      }

//    long d1 = (d_red + 1) / 2;       /* d1 is about m/4 */
//    if (d1 < 1) d1 = 1;
//    if (d1 >= d_red) d1 = d_red - 1;

//    long n = degA - 2 * d1 + 2;      /* n is the ignored part, about m/2 */
//    if (n < 0) n = 0;

//    T A1, B1;

//    if (n != 0)
//      {
//        trunc(A1, A, n);          // A1 = A mod x^n
//        trunc(B1, B, n);          // B1 = B mod x^n
//        RightShift (A, A, n);     // A = A div x^n has about m/2 bits
//        RightShift (B, B, n);     // B = B div x^n has about m/2 bits
//      }

//    TMatrix<T> M1;

//    HALF_GCD (M1, A, B, d1, crossover);
//    // the entries of M1 have m/4 bits, and A,B have been reduced to m/4 bits

//    if (n != 0)
//      {
//        LeftShift (A, A, n);
//        LeftShift (B, B, n);
//        mul (A1, B1, M1);       // A1,B1: m/2  M1:m/4  cost: 4 M(m/2,m/4)
//        add (A, A, A1);
//        add (B, B, B1);
//      }

//    /* now A, B have 3m/4 bits */

//    long d2 = deg(B) - (degA - d_red); /* should be about m/2 */

//    if (IsZero(B) || d2 <= 0)
//      {
//        M_out = M1;
//        return;
//      }

//    T Q;
//    TMatrix<T> M2;

//    DivRem(Q, A, A, B);
//    swap(A, B);

//    T t(INIT_SIZE, deg(M1(1,1))+deg(Q)+1);

//    mul(t, Q, M1(1,0));
//    sub(t, M1(0,0), t);
//    swap(M1(0,0), M1(1,0));
//    swap(M1(1,0), t);

//    t.kill();

//    t.SetMaxLength(deg(M1(1,1))+deg(Q)+1);

//    mul(t, Q, M1(1,1));
//    sub(t, M1(0,1), t);
//    swap (M1(0,1), M1(1,1));
//    swap (M1(1,1), t);

//    if (IsZero(B) || deg(B) <= degA - d_red)
//      {
//        M_out = M1;
//        return;
//      }

//    n = deg(A) - 2 * d2 + 2;     // should be about m/4
//    if (n < 0) n = 0;

//    if (n != 0)
//      {
//        trunc(A1, A, n);
//        trunc(B1, B, n);        // A1,B1 have m/4 bits
//        RightShift(A, A, n);    // A,B have m/2 bits
//        RightShift(B, B, n);
//      }

//    HALF_GCD (M2, A, B, d2, crossover);
//    // now A,B have m/4 bits, like the entries of M2

//    if (n != 0)
//      {
//        LeftShift (A, A, n);
//        LeftShift (B, B, n);
//        mul (A1, B1, M2);
//        add (A, A, A1);
//        add (B, B, B1);
//      }

//    // TODO: use Strassen mul
//    //   mul (M_out, M2, M1); 
//    strassen_mul(M_out, M2, M1);
// }



// -----------------------------------------------------------------------------

// Compute G = GCD(A,B) = U*A + V*B using a subquadratic gcd algorithm.
// The crossover parameter is for profiling purpose only!
template < class T >
void HXGCD(T&G, T& U, T& V, const T& A, const T& B)
{
  T AA(A), BB(B);

  if (IsZero(B)) {
    set(U);
    clear(V);
    G = A;
    return;
  }

  if (IsZero(A)) {
    clear(U);
    set(V);
    G = B;
    return;
  }


  TMatrix<T> M;
  HALF_GCD(M, AA, BB, deg(AA)+1);

  //  G = M(0,0)*A + M(0,1)*B;
  T temp;
  mul(G,M(0,0),A);
  mul(temp,M(0,1),B);
  add(G,G,temp);
  
  T cl;
  cl = inv(LeadCoeff(G));
  mul(U,M(0,0),cl);
  mul(V,M(0,1),cl);
  mul(G,G,cl);
}



//
// HXGCD_LEFT
//

template < class T >
void HXGCD_LEFT(T & G, T & X, const T & A, const T & B)
{
  T Y;
  HXGCD(G,X,Y,A,B);
}



//
// Partial Euclidean algorithm (for NUCOMP)
//

template < class T >
void HXGCD_PARTIAL(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag)
{  
  /*
  TMatrix<T> M_out;

  // if deg(R2) < deg(R1), first quotient is zero.  Swap.
  bool swapped = false;
  if (deg(R2) < deg(R1)) {
    swap(R2,R1);
    swapped = true;
  }

  long d_red = deg(R2) - bound;
  HALF_GCD(M_out, R2, R1, d_red);

  // swap B2 and B1, because of initial swap
  if (swapped) {
    C2 = M_out(1,1);
    C1 = M_out(0,1);
  }
  else {
    C2 = M_out(0,1);
    C1 = M_out(1,1);
  }
  */

  XGCD_PARTIAL_ITER(R2,R1,C2,C1,bound,flag);
}



//
// Partial Euclidean algorithm (for fast reduce) with revised termination for 
// polynomial types
//

template < class T >
void HXGCD_PARTIAL_REDUCE(T & R2, T & R1, T & C2, T & C1, long bound, bool & flag, bool even)
{  
  /*
  TMatrix<T> M_out;

  // if deg(R2) < deg(R1), first quotient is zero.  Swap.
  bool swapped = false;
  if (deg(R2) < deg(R1)) {
    swap(R2,R1);
    swapped = true;
  }

  long d_red = deg(R2) - bound;

  HALF_GCD(M_out, R2, R1, d_red);


  // swap B2 and B1, because of initial swap
  if (swapped) {
    B2 = M_out(0,1);
    B1 = M_out(1,1);
  }
  else {
    B2 = M_out(0,0);
    B1 = M_out(1,0);
  }

  if (deg(R1) == bound && !even) {
    T Q,t;

    DivRem(Q, R2, R2, R1);
    swap(R1, R2);

      // r = B2 + q B1
      mul(t,Q,B1);
      add(t,B2,t);
      B2 = B1;
      B1 = t;
  }
  */

  XGCD_PARTIAL_REDUCE_ITER(R2,R1,C2,C1,bound,flag,even);
}
