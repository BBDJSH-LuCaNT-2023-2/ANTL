#ifndef ANTL_VECTOR_FUNCTIONS_HPP
#define ANTL_VECTOR_FUNCTIONS_HPP

#include "../Arithmetic/QQ.hpp"
#include "../common.hpp"
#include <ANTL/Cubic/CubicIdeal.hpp>
#include <vector>
using namespace NTL;
using namespace ANTL;


#include <vector>
/**
* @brief
* @param[in]
* @param[out]
*/
  template<typename PType>
  bool positive_func(PType dubby);
  /**
  * @brief Checks whether all coordinates of vec1 and vec2 are within maxdist of each other
  * Currently for RR and doubles
  * @param[out] true if abs(v1[i] -v2[i]) < maxdist for all i
  * @param[in] B vec1 is a vector of real type entries.
  * @param[in] vec2 is a vector with real type entries
  * @param[in] maxdist is a real type value that indicates the max acceptable distance of entries to be considered close
  * @pre vec1 and vec2 need to be the same size.
  */


  template<typename PType>
  bool is_close(const std::vector<PType> & vec1, const std::vector<PType> & vec2, const PType & maxdist);


  template<typename Type, typename PType>
  void compute_initial_s(const std::vector<PType> & alpha, const int kbound);

  /*
  * @brief Takes a length r log vector and returns the corresponding length r+1 exponentiated vector. equivalent to create_target from the pari stuff.
  * @param[in] log_vec is a vector with real entries having length equal to the unit rank r of the field
  * @param[out] valuationvec is a length r+1 vector, consisting of  exp(log_vec[i]), -sum(deg(i)*log_vec[i])/deg(r+1) )
  * Here deg(i) is 1 if the ith coordinate is a real embedding, 2 otherwise.
  */
  template<typename PType>
  void log_to_valuation(std::vector<PType> &valuationvec, const std::vector<PType> & log_vec, const int r1);

#include "VectorFunctions_impl.hpp"
#endif
