/**
 * @file Character.hpp
 * @author Leonard Nooy, Mike Jacobson
 * @version $Header$
 *
 * This file contains the class definition for a simple character class.
 * It includes functions for quadartic and quartic characters.
 */

#ifndef ANTL_L_FUNCTION_CHARACTER_H
#define ANTL_L_FUNCTION_CHARACTER_H

#include <ANTL/common.hpp>
#include <ANTL/Arithmetic/CC.hpp>

// Forward Declarations
namespace NTL
{
  class GF2EX;
  class GF2EXModulus;
  class ZZ_pX;
  class ZZ_pEX;
  class zz_pX;
  class zz_pEX;
}

namespace ANTL
{

#define QUADRATIC_MODE 2
#define QUARTIC_MODE   4


  /**
   * definition of Character class (character computation using precomputed
   * values)
   */

  template < class T >
  class Character
  {
  protected:
    // used by all characters
    T M;                          // character modulus
    T h;                          // h required for GF2EX (y^2 + hy = f model)
    int mode;                     // initialization mode:
    //   0 - uninitialized
    //   2 - quadratic
    //   4 - quartic

    // used by quadratic characters
    int D2;                       // M mod 2
    int D4;                       // M mod 4
    int D8;                       // M mod 8

    // used by quartic characters
    long rootp;			// sqrt(p)
    T MPList[4];			// this holds the precalculated values of 
    //   Mp = 2^(p-1)/4 mod(p)
    CC < float >CHI1List[4];	// holds pre-calculated values of chi(-1)
    CC < float >IList[4];		// holds the values of i^x, where x = {0,1,2,3}
    CC < float >*CHIList;		// holds all values of chi <= sqrt(p);
    long CHIListSize;		// the number of elements of CHIList


  public:
    /* Constructor/Destructor */
    Character ();
    ~Character ();

    /* Initalization function */
    bool init(const T & newM, int new_mode);
    bool init(const T & newM, const T & new_h, int new_mode);

    /* Quadratic character functions */
    long quadratic(const T & n);
    long quadratic_long(const long & n);       // to be implemented for T=longlong and ZZ

    /* Quartic character functions */
    bool allocate_quartic_table();
    bool init_quartic_table ();
    CC < float > quartic(const T & n);
    CC < float > quartic_table(const T & n);
    CC < float > quartic_long(const long & n);
    CC < float > quartic_table_long(const long & n);
  };

  //
  // Specialized Method Prototypes
  //

  // initialization
  template <> bool Character<long>::init(const long & newM, int new_mode);
  template <> bool Character<long long>::init(const long long & newM, int new_mode);
  template <> bool Character<ZZ>::init(const ZZ & newM, int new_mode);
  template <> bool Character<GF2EX>::init(const GF2EX & newM, const GF2EX & newh, int new_mode);

  // quadratic
  template <> long Character<GF2EX>::quadratic(const GF2EX & n);
  template <> long Character<long long>::quadratic_long(const long & n);
  template <> long Character<ZZ>::quadratic_long(const long & n);

  // quartic
  template <> bool Character<ZZ>::allocate_quartic_table ();
  template <> bool Character<long>::allocate_quartic_table ();
  template <> bool Character<long long>::allocate_quartic_table ();
  template <> bool Character<ZZ>::init_quartic_table ();
  template <> bool Character<long>::init_quartic_table();
  template <> bool Character<long long>::init_quartic_table();
  template <> CC < float > Character<ZZ>::quartic (const ZZ & n);
  template <> CC < float > Character<long>::quartic (const long & n);
  template <> CC < float > Character<long long>::quartic (const long long & n);
  template <> CC < float > Character<ZZ>::quartic_long (const long & n);
  template <> CC < float > Character<long long>::quartic_long (const long & n);
  template <> CC < float > Character<ZZ>::quartic_table (const ZZ & n);
  template <> CC < float > Character<long>::quartic_table (const long & n);
  template <> CC < float > Character<long long>::quartic_table (const long long & n);
  template <> CC < float > Character<ZZ>::quartic_table_long (const long & n);
  template <> CC < float > Character<long long>::quartic_table_long (const long & n);


} // ANTL

#include "../../../src/L_function/Character_impl.hpp" 

#endif // guard
