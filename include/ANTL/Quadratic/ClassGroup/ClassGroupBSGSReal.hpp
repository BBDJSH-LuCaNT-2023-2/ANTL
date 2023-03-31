#ifndef CLASSGROUP_BSGS_H
#define CLASSGROUP_BSGS_H

#include <ANTL/common.hpp>
#include <unordered_set>

#include <ANTL/HashTable/HashEntryInt.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>

#include <ANTL/LinearAlgebra/Smith.hpp>

#include <ANTL/Quadratic/QuadraticClassGroupElement.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

// template <> struct std::hash<QuadraticIdealBase<ZZ>> {
//   std::size_t operator()(QuadraticIdealBase<ZZ> const &qib) const noexcept {
//     std::size_t h1 = std::hash<int>{}(to<int>(qib.get_a()));
//     std::size_t h2 = std::hash<int>{}(to<int>(qib.get_b()));
//     return h1 ^ (h2 << 1); // or use boost::hash_combine
//   }
// };

NTL_CLIENT;
using namespace ANTL;

namespace ANTL {

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class T> class ClassGroupBSGSReal {
private:

  // DBG_CONSTANTS
  bool DBG_CGBGRL = false;
  bool DBG_ISPRIN = false;
  bool DBG_GNEXTP = false;

  QuadraticOrder<T> *quadratic_order;

  double regulator;

  T delta;

  // is_principal variables
  double sqrt_regulator;
  std::unordered_set<QuadraticIdealBase<T>> baby_step_list;
  QuadraticClassGroupElement<T> is_principal_giant_step;

  // factor base for computing CL
  QuadraticClassGroupElement<T> *fact_base;

  // indice of fact base elements that
  long *contributors;

  // number of primes used for BSGS
  long num_prime_ideals;

  // class number
  ZZ h;

  // invariants of CL
  vector<ZZ> CL;

  // size of fact_base
  long numFB;

  // transformation matrix
  mat_ZZ U_mat;

  // A Toggle for resetting a PrimeSeq.
  // This seems to be neccessary since sometimes assign_prime(G) can return not
  // just a principal ideal but the identity itself; which on the next call of
  // get_next_prime() (that is, get_next_prime(assign_prime(G))) can reset PS
  // and cause an infinite loop. For the exact code being referred to, see the
  // do-while loop immediately below CGBGRL: STEP 4.1
  bool reset_prime_seq;

public:
  ClassGroupBSGSReal(QuadraticOrder<T> *quadratic_order);

  ~ClassGroupBSGSReal();

  void cg_bsgs_real(const ZZ &hstar);

  void set_regulator(double &ext_regulator){
    regulator = ext_regulator;
    sqrt_regulator = sqrt(regulator);
  }

  vector<ZZ> get_class_group(){return CL;}

private:
  ZZ get_dist_mod(const T &Delta) { return CeilToZZ(log(to_RR(to_ZZ(Delta)))); }

  void get_next_prime(QuadraticClassGroupElement<T> &G);

  void is_principal_init();

  bool is_principal(const QuadraticClassGroupElement<T> &G);

  void decode_vector(mat_ZZ &Bmat, const ZZ &Bjj, const ZZ &r, const ZZ &q,
                    vec_ZZ &Rvec, vec_ZZ &Qvec, long nR, long nQ);
};

} // ANTL

#endif
