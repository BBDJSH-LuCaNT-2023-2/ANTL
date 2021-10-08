#ifndef ANTL_VECTOR_FUNCTIONS_HPP
#define ANTL_VECTOR_FUNCTIONS_HPP

using namespace ANTL;
using namespace NTL;

#include <vector>


/**
* @brief Check if two vectors have all entries close together. Currently for RR and doubles
* @param[out] true if abs(v1[i] -v2[i]) < maxdist for all i
* @param[in] B vec1 is a vector of real type entries.
* @param[in] vec2 is a vector with real type entries
* @param[in] maxdist is a real type value that indicates the max acceptable distance of entries to be considered close
* @pre vec1 and vec2 need to be the same size.
*/
template<typename PType>
bool is_close(std::vector<PType> & vec1, std::vector<PType> & vec2, PType & maxdist);

#endif
