/**
 * @file Character_impl.hpp
 * @author Leonard Nooy / Mike Jacobson
 * @remarks This file is to be included from Character.h only.
 * @version $Header$
 *
 * This file contains the code that the compiler will instantiate when the
 * user attempts to create a class that is not one of the types: long,
 * long long, ZZ, or GF2EX
 */

namespace ANTL
{

  /**
   * Implementations of Character class member functions
   */

  /*
   * Function: Character::Character
   * Description: Contructor
   * Inputs: void
   * Outputs: void
   */
  template < class T >
  Character<T>::Character ()
  {
    ::clear(M);
    mode = 0;

    /* initalize the IList */
    IList[0].set ();
    IList[1].set_i ();
    IList[2].negate (IList[0]);
    conjugate(IList[3], IList[1]);

    CHIList = 0;
    CHIListSize = 0;
  }


  /*
   * Function: Character::~Character
   * Description: Destructor
   * Inputs: void
   * Outputs: void
   */
  template < class T >
  Character<T>::~Character ()
  {
    if (CHIList)
      {
	delete[]CHIList;
	CHIList = 0;
      }
  }



  /*
   * Function: Character::init(M,mode)
   * Description: This function must be called before calling and of the 
   *   character functions
   * Inputs: T newM - new character modulus
   *         int mode - mode (2 - quadratic, 4 - quartic)
   * Outputs: bool - true = success
   *               - false = failure
   */
  template < class T >
  bool Character<T>::init(const T & newM, int new_mode)
  {
    M = newM;
    mode = new_mode;

    if (mode == QUARTIC_MODE) {
      cerr << "Quartic characters undefined for poly base types" << endl;
      mode = 0;
      return false;
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
  template < class T >
  bool Character<T>::init(const T & newM, const T & newh, int new_mode)
  {
    cerr << "Use init(M,mode) for this base type." << endl;
    mode = 0;
    return false;
  }



  //
  // quadratic characters
  //

  template < class T >
  long Character<T>::quadratic(const T & n)
  {
    return Kronecker(M,n);
  }



  template < class T >
  long Character<T>::quadratic_long(const long & n)
  {
    return 0;
  }



  //
  // quartic characters
  //

  template < class T >
  bool Character<T>::allocate_quartic_table()
  {
    return false;
  }



  /*
   * Function::init_quartic_table
   * Description: This function must be called before using the quartic_table function
   * Inputs: void.
   * Outputs: bool - true = success
   *                 false = failure
   */

  template < class T >
  bool Character<T>::init_quartic_table()
  {
    return false;
  }



  /*
   * Function::quartic
   * Description: This function calculates the value of the imaginary quartic
   *   character using the technique presented by Stephen Louboutin, but using
   *   some pre-calculated tables.  Assumes the Character class has been properly
   *   initialized using init(M,mode).
   * Inputs: ZZ a - The value of (a/D)_4
   * Outputs: CC<float> - The value of the character
   */
  template < class T >
  CC < float >
  Character<T>::quartic(const T & n)
  {
    CC < float > x;
    x.clear();
    return x;
  }

  template < class T >
  CC < float >
  Character<T>::quartic_long(const long & n)
  {
    CC < float > x;
    x.clear();
    return x;
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
  template < class T >
  CC < float >
  Character<T>::quartic_table(const T & n)
  {
    CC <float> x;
    x.clear();
    return x;
  }

  template < class T >
  CC < float >
  Character<T>::quartic_table_long(const long & n)
  {
    CC <float> x;
    x.clear();
    return x;
  }

} // ANTL

