#ifndef QUADINDCALC_IMPL_H
#define QUADINDCALC_IMPL_H

// this file has stub implementations for compute_cg_fac_base, compute_cg_relations, compute_cg_mat
// Calling these functions in the right order is handled by the parent class
// Computations are passed off to the classes FactorBase and RelationGenerator

#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/QuadRelationGenerator.hpp"
#include "ANTL/IndexCalculus/Relation/QuadRelation.hpp"
#include "ANTL/IndexCalculus/FactorBase/QuadFactorBase.hpp"

template <class T, class R>
void QuadIndCalc<T,R>::compute_cg_fac_base(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  //cf ANTL/include/ANTL/quadratic/index_calculus/factor_base.hpp
  std::string fb_size_str = "size_fb";
  if ( params.find(fb_size_str) == params.end() ) {
    NTL::ZZ fb_size(NTL::INIT_VAL, params.find(fb_size_str)->first);
    IndCalc<R>::cg_factor_base = QuadFactorBase<T>(fb_size);

    // get factor basis elements using QuadFactorBase
  } else {
    std::cout << "compute_cg_fac_base: params should specify factor basis size as fb_size" << std::endl;
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_cg_relations(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  //cf ANTL/src/quadratic/qo_relation_generator_impl.hpp
  std::string num_relations_str = "num_relations";
  if ( params.find(num_relations_str) == params.end() ) {
    NTL::ZZ num_relations(NTL::INIT_VAL, params.find(num_relations_str)->first);
    IndCalc<R>::cg_relation_generator = QuadRelationGenerator<T>(quad_order, params, IndCalc<R>::cg_factor_base);

    // find relations by calling RelationGenerator's get_relation function
    // Call this function num_relations times

  } else {
    std::cout << "compute_cg_relations: params should specify factor basis size as fb_size" << std::endl;
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_cg_mat(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  //cf ANTL/include/ANTL/quadratic/index_calculus/qo_relation_matrix.hpp
  IndCalc<R>::cg_mat.SetDims(IndCalc<R>::cg_relations.size(), IndCalc<R>::cg_factor_base.size());

  // set matrix using functions in QuadRelation
};

#endif //QUADCLASSGROUPINDCALC_H
