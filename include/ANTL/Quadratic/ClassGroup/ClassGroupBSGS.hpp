#ifndef CLASSGROUP_BSGS_H
#define CLASSGROUP_BSGS_H

#include <ANTL/common.hpp>

#include <ANTL/HashTable/HashEntryReal.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>
#include <ANTL/L_function/L_function.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>
#include <ANTL/Quadratic/QuadraticOrder.hpp>

NTL_CLIENT;
using namespace ANTL;

//DBG_CONSTANTS
bool DBG_CGBSGS = false;

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class T, class U> class ClassGroupBSGS {
private:

public:

private:

};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.


#endif
