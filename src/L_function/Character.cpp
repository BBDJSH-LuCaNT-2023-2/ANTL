/**
 * @file Character.cpp
 * @version $Header$
 */

#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>

#include <ANTL/L_function/Character.hpp>

using namespace std;

namespace ANTL
{

  /**
   * Implementations of Character class member functions
   */

  /*
   * Function: Character::init(M,mode)
   * Description: This function must be called before calling and of the 
   *   character functions
   * Inputs: T newM - new character modulus
   *         int mode - mode (2 - quadratic, 4 - quartic)
   * Outputs: bool - true = success
   *               - false = failure
   */
  template <>
  bool Character<long>::init(const long & newM, int new_mode)
  {
    M = newM;
    mode = new_mode;

    if (mode == QUADRATIC_MODE) {
      D2 = M & 1;
      D4 = M & 3;
      D8 = M & 7;
    }
    else if (mode == QUARTIC_MODE) {
      /* calculate MPList */
      long Mp;

      D8 = M & 7;

      // make sure that prime is congruent to 5 mod 8
      if (D8 != 5)
	{
	  cout << "(p % 8) != 5, please use a valid prime." << endl;
	  return false;
	}

      Mp = PowerMod (2, (M - 1) >> 2, M);
      MPList[0] = 1;
      MPList[1] = Mp;
      MPList[2] = MulMod(Mp,Mp,M);
      MPList[3] = MulMod(MPList[2],Mp,M);
    }

    return true;
  }


  template <>
  bool Character<long long>::init(const long long & newM, int new_mode)
  {
    M = newM;
    mode = new_mode;

    if (mode == QUADRATIC_MODE) {
      D2 = M & 1;
      D4 = M & 3;
      D8 = M & 7;
    }
    else if (mode == QUARTIC_MODE) {
      cerr << "Quartic characters undefined for long long  base types" << endl;
      mode = 0;
      return false;

      /* calculate the precompiled MPList for usage by chi1 */
      long Mp;

      // make sure that prime is congruent to 5 mod 8
      if ((M & 7) != 5)
	{
	  cout << "(p % 8) != 5, please use a valid prime." << endl;
	  return false;
	}

      Mp = PowerMod ((long long) 2, (M - 1) >> 2, M);
      MPList[0] = 1;
      MPList[1] = Mp;
      MPList[2] = MulMod(Mp,Mp,M);
      MPList[3] = MulMod(MPList[2],Mp,M);
    }

    return true;
  }


  template <>
  bool Character<ZZ>::init(const ZZ & newM, int new_mode)
  {
    M = newM;
    mode = new_mode;

    if (mode == QUADRATIC_MODE) {
      D2 = IsOdd(M);
      D4 = trunc_long(M,2);
      D8 = trunc_long(M,3);
      if (M < 0) {
	D4 = (4 - D4);
	D8 = (8 - D8);
      }
    }
    else if (mode == QUARTIC_MODE) {
      /* calculate the precompiled MPList for usage by chi1 */
      ZZ Mp;				// temp variable

      D8 = trunc_long(M,3);
      if (M < 0)
	D8 = (8 - D8);

      // make sure that prime is congruent to 5 mod 8
      if (D8 != 5)
	{
	  cout << "(p % 8) != 5, please use a valid prime." << endl;
	  return false;
	}

      PowerMod (Mp, to_ZZ (2), (M - 1) >> 2, M);
      MPList[0] = 1;
      MPList[1] = Mp;
      SqrMod(MPList[2], Mp, M);
      MulMod(MPList[3], MPList[2], Mp, M);
    }

    return true;
  }


  /*
   * Function: Character::init(M,h,mode)
   * Description: This function must be called before calling and of the 
   *   character functions
   * Inputs: T newM - new character modulus
   *         int mode - mode (2 - quadratic, 4 - quartic)
   * Outputs: bool - true = success
   *               - false = failure
   */
  template <>
  bool Character<GF2EX>::init(const GF2EX & newM, const GF2EX & newh, int new_mode)
  {
    M = newM;
    h = newh;
    mode = new_mode;

    if (mode == QUARTIC_MODE) {
      cerr << "Quartic characters undefined for poly base types" << endl;
      mode = 0;
      return false;
    }

    return true;
  }



  //
  // quadratic characters
  //

  template <>
  long Character<long>::quadratic(const long & n)
  {
    long CHI;

    if (!(n & 1) && !D2) // a and n are both even
      CHI = 0;
    else
      {
	bool neg = false;
	long t = n;
	long temp;

	while (!(t & 1))
	  {
	    t >>= 1;
	    neg = !neg;
	  }

	temp = M % t;
	if (temp < 0)  temp += t;
	CHI = Jacobi_base(temp,t);
	if (neg && (D8 == 5))
	  CHI = -CHI;
      }

    return CHI;
  }


  template <>
  long Character<long long>::quadratic(const long long & n)
  {
    long CHI;

    if (!(n & 1) && !D2) // D and p are both even
      CHI = 0;
    else
      {
	bool neg = false;
	long long t = n;
	long long temp;

	while (!(t & 1))
	  {
	    t >>= 1;
	    neg = !neg;
	  }

	temp = M % t;
	if (temp < 0)  temp += t;
	CHI = Jacobi_base(temp,t);
	if (neg && (D8 == 5))
	  CHI = -CHI;
      }

    return CHI;
  }


  template <>
  long Character<ZZ>::quadratic(const ZZ & n)
  {
    long CHI;

    if (!IsOdd(n) && !D2) // D and p are both even
      CHI = 0;
    else
      {
	bool neg = false;
	ZZ t = n;
	ZZ temp;

	while (!IsOdd(t))
	  {
	    t >>= 1;
	    neg = !neg;
	  }

	rem(temp,M,t);
	CHI = Jacobi_base(temp,t);
	if (neg && (D8 == 5))
	  CHI = -CHI;
      }

    return CHI;
  }


//   template <>
//   long Character<GF2EX>::quadratic(const GF2EX & n)
//   {
//     return Jacobi(h,M,n);
//   }



  template <>
  long Character<long long>::quadratic_long(const long & n)
  {
    long CHI;

    if (!(n & 1) && !D2) // D and p are both even
      CHI = 0;
    else
      {
	bool neg = false;
	long t = n;
	long temp;

	while (!(t & 1))
	  {
	    t >>= 1;
	    neg = !neg;
	  }

	temp = M % t;
	if (temp < 0)  temp += t;
	CHI = Jacobi_base(temp,t);
	if (neg && (D8 == 5))
	  CHI = -CHI;
      }

    return CHI;
  }


  template <>
  long Character<ZZ>::quadratic_long(const long & n)
  {
    long CHI;

    if (!(n & 1) && !D2) // a and n are both even
      CHI = 0;
    else
      {
	bool neg = false;
	long t = n;
	long temp;

	while (!(t & 1))
	  {
	    t >>= 1;
	    neg = !neg;
	  }

	temp = rem(M,t);
	CHI = Jacobi_base(temp,t);
	if (neg && (D8 == 5))
	  CHI = -CHI;
      }

    return CHI;
  }



  //
  // quartic characters
  //

  template <>
  bool Character<ZZ>::allocate_quartic_table ()
  {
    if (CHIListSize < SqrRoot (M) || !CHIList)
      {
	if (CHIList)
	  {
	    delete[]CHIList;
	    CHIList = 0;
	  }

	conv (CHIListSize, SqrRoot (M));
	CHIList = new CC < float >[CHIListSize + 1];

	if (!CHIList)
	  {
	    cout << "Error Allocating CHIList, Characters.h";
	    return false;
	  }
      }

    return true;
  }

  template <>
  bool Character<long>::allocate_quartic_table ()
  {
    if (CHIListSize < SqrRoot(M) || !CHIList)
      {
	if (CHIList)
	  {
	    delete[]CHIList;
	    CHIList = 0;
	  }

	CHIListSize = SqrRoot(M);
	CHIList = new CC < float >[CHIListSize + 1];

	if (!CHIList)
	  {
	    cout << "Error Allocating CHIList, Characters.h";
	    return false;
	  }
      }

    return true;
  }

//   template <>
//   bool Character<long long>::allocate_quartic_table ()
//   {
//     if (CHIListSize < SqrRoot (M) || !CHIList)
//       {
// 	if (CHIList)
// 	  {
// 	    delete[]CHIList;
// 	    CHIList = 0;
// 	  }
//
// 	CHIListSize = SqrRoot (M);
// 	CHIList = new CC < float >[CHIListSize + 1];
//
// 	if (!CHIList)
// 	  {
// 	    cout << "Error Allocating CHIList, Characters.h";
// 	    return false;
// 	  }
//       }
//
//     return true;
//   }


  /*
   * Function::init_quartic_table
   * Description: This function must be called before using the quartic_table function
   * Inputs: void.
   * Outputs: bool - true = success
   *                 false = failure
   */

  template <>
  bool Character<ZZ>::init_quartic_table ()
  {
    long i;

    // make sure that prime is congruent to 5 mod 8
    if (D8 != 5 || (mode != QUARTIC_MODE))
      {
	cout << "Please call init with a valid prime and quartic mode before calling this function" << endl;
	return false;
      }

    /* Calculate the values for (-1 / p)_4 */
    CHI1List[1] = quartic(M - 1);
    CHI1List[0].set ();
    square(CHI1List[2], CHI1List[1]);
    multiply(CHI1List[3], CHI1List[2], CHI1List[1]);

    /* Now precompute a values of Chi from 1 -> sqrt(p) for chi3 */
    if (!allocate_quartic_table())
      return false;

    CHIList[0].clear ();

    rootp = to_long (SqrRoot (M));
    i = 1;
    while (i <= rootp)
      {
	CHIList[i] = quartic_long (i);
	i++;
      }

    return true;
  }


  template <>
  bool Character<long>::init_quartic_table()
  {
    long i;

    // make sure that prime is congruent to 5 mod 8
    if ((M & 7) != 5 || (mode != QUARTIC_MODE))
      {
	cout << "Please call init with a valid prime and quartic mode before calling this function" << endl;
	return false;
      }

    /* Calculate the values for (-1 / p)_4 */
    CHI1List[1] = quartic (M - 1);
    CHI1List[0].set ();
    square(CHI1List[2], CHI1List[1]);
    multiply(CHI1List[3], CHI1List[2], CHI1List[1]);

    /* Now precompute a values of Chi from 1 -> sqrt(p) for chi3 */
    if (!allocate_quartic_table())
      return false;

    CHIList[0].clear ();

    rootp = SqrRoot (M);
    i = 1;
    while (i <= rootp)
      {
	CHIList[i] = quartic (i);
	i++;
      }

    return true;
  }


//   template <>
//   bool Character<long long>::init_quartic_table()
//   {
//     long i;
//
//     make sure that prime is congruent to 5 mod 8
//     if ((M & 7) != 5 || (mode != QUARTIC_MODE))
//       {
// 	cout << "Please call init with a valid prime and quartic mode before calling this function" << endl;
// 	return false;
//       }
//
//     /* Calculate the values for (-1 / p)_4 */
//     CHI1List[1] = quartic (M - 1);
//     CHI1List[0].set ();
//     square(CHI1List[2], CHI1List[1]);
//     multiply(CHI1List[3], CHI1List[2], CHI1List[1]);
//
//     /* Now precompute a values of Chi from 1 -> sqrt(p) for chi3 */
//     if (!allocate_quartic_table())
//       return false;
//
//     CHIList[0].clear ();
//
//     rootp = SqrRoot (M);
//     i = 1;
//     while (i <= rootp)
//       {
// 	CHIList[i] = quartic (i);
// 	i++;
//       }
//
//     return true;
//   }



  /*
   * Function::quartic
   * Description: This function calculates the value of the imaginary quartic
   *   character using the technique presented by Stephen Louboutin, but using
   *   some pre-calculated tables.  Assumes the Character class has been properly
   *   initialized using init(M,mode).
   * Inputs: ZZ a - The value of (a/D)_4
   * Outputs: CC<float> - The value of the character
   */
  template <>
  CC < float >
  Character<ZZ>::quartic (const ZZ & n)
  {
    register long i;
    ZZ n4;
    ZZ temp;

    rem(temp,n,M);
    if (IsZero(temp))
      {
	CC < float >x;

	x.clear ();
	return x;
      }

    PowerMod (n4,temp, (M - 1) >> 2, M);

    for (i=0; MPList[i] != n4; ++i);
    return IList[i];
  }

  template<>
  CC < float >
  Character<long>::quartic (const long & n)
  {
    register long i;
    long n4;
    long temp;

    temp = n % M;
    if (temp == 0)
      {
	CC < float >x;

	x.clear ();
	return x;
      }

    n4 = PowerMod (temp, (M - 1) >> 2, M);

    for (i=0; MPList[i] != n4; ++i);
    return IList[i];
  }

//   template<>
//   CC < float >
//   Character<long long>::quartic (const long long & n)
//   {
//     register long i;
//     long long n4;
//     long long temp;
//
//     temp = n % M;
//     if (temp == 0)
//       {
// 	CC < float >x;
//
// 	x.clear ();
// 	return x;
//       }
//
//     n4 = PowerMod (temp, (M - 1) >> 2, M);
//
//     for (i=0; MPList[i] != n4; ++i);
//     return IList[i];
//   }


  template <>
  CC < float >
  Character<ZZ>::quartic_long (const long & n)
  {
    register long i;
    ZZ n4;
    ZZ temp;

    rem(temp,to_ZZ(n),M);
    if (IsZero(temp))
      {
	CC < float >x;

	x.clear ();
	return x;
      }

    PowerMod (n4,temp, (M - 1) >> 2, M);

    for (i=0; MPList[i] != n4; ++i);
    return IList[i];
  }

  template <>
  CC < float >
  Character<long>::quartic_long (const long & n)
  {
    register long i;
    long long n4;
    long long temp;

    temp = n % M;
    if (!temp)
      {
	CC < float >x;

	x.clear ();
	return x;
      }

    n4 = PowerMod (temp, (M - 1) >> 2, M);

    for (i=0; MPList[i] != n4; ++i);
    return IList[i];
  }



  /*
   * Function::quartic_table
   * Description: This function calculates the value of the imaginary quartic 
   *  character using the technique presented by Micheal Jacobson. It uses the
   *  concept of the euclidien algorithm to find two smaller factors who's
   *  values are looked up in a table.  Assumes character class has been 
   *  properly initialized using init(M,mode) and that the table has been
   *  initialized using init_quartic_table().
   * Inputs: ZZ a - The value of (a/D)_4
   * Outputs: CC<float> - The value of the character
   */
  template <>
  CC < float >
  Character<ZZ>::quartic_table (const ZZ & n)
  {
    CC < float >x, y;
    ZZ B_2, B_1, B;
    ZZ R_2, R_1, R;
    ZZ q;
    long i;

    // if n <= sqrt(p) then just look up the value of Chi(a) from the table
    if (n <= rootp)
      {
	return CHIList[to_long (n)];
      }

    /*
      PRECONDITION:   a not divisible by p!

      else if((a % p) == 0){
      x.clear();
      return x;
      }
    */

    i = 0;

    R_2 = n;
    DivRem(B,R_1,M,n);

    set(B_2);
    B_1 = B;

    while (B <= rootp)
      {
	// use the euclidean algorithm to find R, B such that X(a) = (-1/p)^i(R/p)(B/p)^-1
	DivRem (q, R, R_2, R_1);

	R_2 = R_1;
	R_1 = R;

	mul(B,B_1,q);
	add(B,B,B_2);
	B_2 = B_1;
	B_1 = B;

	i++;
      }

    conjugate(y, CHIList[to_long(B_2)]);
    multiply(x, y, CHIList[to_long (R_2)]);
    multiply(x, x, CHI1List[i & 3]);
    return x;
  }

  template <>
  CC < float >
  Character<long>::quartic_table (const long & n)
  {
    CC < float >x, y;
    long B_2, B_1, B;
    long R_2, R_1, R;
    long q;
    long i;

    // if a <= sqrt(p) then just look up the value of Chi(a) from the table
    if (n <= rootp)
      {
	return CHIList[n];
      }

    /*
      PRECONDITION:   a nod divisible by p!

      else if((a % plong) == 0){
      x.clear();
      return x;
      }
    */

    i = 0;

    R_2 = n;
    DivRem(B,R_1,M,n);

    B_2 = 1;
    B_1 = B;

    while (B <= rootp)
      {
	// use the euclidean algorithm to find R, B such that X(a) = (-1/p)^i(R/p)(B/p)^-1
	DivRem(q,R,R_2,R_1);
	R_2 = R_1;
	R_1 = R;

	B = (q * B_1) + B_2;
	B_2 = B_1;
	B_1 = B;

	i++;
      }

    conjugate(y, CHIList[B_2]);
    multiply(x, y, CHIList[R_2]);
    multiply(x, x, CHI1List[i & 3]);
    return x;
  }


//   template <>
//   CC < float >
//   Character<long long>::quartic_table (const long long & n)
//   {
//     CC < float >x, y;
//     long long B_2, B_1, B;
//     long long R_2, R_1, R;
//     long long q;
//     long i;
//
//     if a <= sqrt(p) then just look up the value of Chi(a) from the table
//     if (n <= rootp)
//       {
// 	return CHIList[n];
//       }
//
//     /*
//       PRECONDITION:   a nod divisible by p!
//
//       else if((a % plong) == 0){
//       x.clear();
//       return x;
//       }
//     */
//
//     i = 0;
//
//     R_2 = n;
//     DivRem(B,R_1,M,n);
//
//     B_2 = 1;
//     B_1 = B;
//
//     while (B <= rootp)
//       {
// 	use the euclidean algorithm to find R, B such that X(a) = (-1/p)^i(R/p)(B/p)^-1
// 	DivRem(q,R,R_2,R_1);
// 	R_2 = R_1;
// 	R_1 = R;
//
// 	B = (q * B_1) + B_2;
// 	B_2 = B_1;
// 	B_1 = B;
//
// 	i++;
//       }
//
//     conjugate(y, CHIList[B_2]);
//     multiply(x, y, CHIList[R_2]);
//     multiply(x, x, CHI1List[i & 3]);
//     return x;
//   }


  template <>
  CC < float >
  Character<ZZ>::quartic_table_long (const long & n)
  {
    CC < float > x, y;
    ZZ B_2, B_1, B;
    long R_2, R_1, R;
    long q;
    long i;

    // if n <= sqrt(p) then just look up the value of Chi(a) from the table
    if (n <= rootp)
      {
	return CHIList[n];
      }

    /*
      PRECONDITION:   a not divisible by p!

      else if((a % p) == 0){
      x.clear();
      return x;
      }
    */

    i = 0;

    R_2 = n;
    R_1 = DivRem(B,M,n);

    set(B_2);
    B_1 = B;

    while (B <= rootp)
      {
	// use the euclidean algorithm to find R, B such that X(a) = (-1/p)^i(R/p)(B/p)^-1
	DivRem (q, R, R_2, R_1);

	R_2 = R_1;
	R_1 = R;

	mul(B,B_1,q);
	add(B,B,B_2);
	B_2 = B_1;
	B_1 = B;

	i++;
      }

    conjugate(y, CHIList[to_long(B_2)]);
    multiply(x, y, CHIList[R_2]);
    multiply(x, x, CHI1List[i & 3]);

    return x;
  }

//   template <>
//   CC < float >
//   Character<long long>::quartic_table_long (const long & n)
//   {
//     CC < float > x, y;
//     long long B_2, B_1, B;
//     long R_2, R_1, R;
//     long q;
//     long i;
//
//     if n <= sqrt(p) then just look up the value of Chi(a) from the table
//     if (n <= rootp)
//       {
// 	return CHIList[n];
//       }
//
//     /*
//       PRECONDITION:   a not divisible by p!
//
//       else if((a % p) == 0){
//       x.clear();
//       return x;
//       }
//     */
//
//     i = 0;
//
//     R_2 = n;
//     R_1 = DivRem(B,M,n);
//
//     set(B_2);
//     B_1 = B;
//
//     while (B <= rootp)
//       {
// 	use the euclidean algorithm to find R, B such that X(a) = (-1/p)^i(R/p)(B/p)^-1
// 	DivRem (q, R, R_2, R_1);
//
// 	R_2 = R_1;
// 	R_1 = R;
//
// 	B = B_1*q + B_2;
// 	B_2 = B_1;
// 	B_1 = B;
//
// 	i++;
//       }
//
//     conjugate(y, CHIList[to<long>(B_2)]);
//     multiply(x, y, CHIList[R_2]);
//     multiply(x, x, CHI1List[i & 3]);
//
//     return x;
//   }

} // ANTL

