/**
 * @file smith.hpp
 * @author Kris Luttmer
 * @version $Header$
 */

#ifndef ANTL_SMITH_H
#define ANTL_SMITH_H

#include <ANTL/common.hpp>
#include <NTL/ZZ.h>
#include <NTL/mat_ZZ.h>

// Forward declarations.
/*
namespace NTL
{
  class mat_ZZ;
}
*/

NTL_CLIENT;
using namespace ANTL;

namespace ANTL
{
  /*****************************************************************************
	An algorithm to convert a matrix to Smith Normal Form based on algorithm 2.4.14
	from "A Course in Computational Algebraic Number Theory" by Henri Cohen
	Written by: Kris Luttmer
	Date: May 14, 2003
  *****************************************************************************/
  mat_ZZ SmithNormalForm (mat_ZZ A, mat_ZZ & U, mat_ZZ & V);

} //nf

#endif // guard
