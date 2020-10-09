#ifndef QUADINDCALC_IMPL_H
#define QUADINDCALC_IMPL_H

// this file has stub implementations for compute_cg_fac_base, compute_cg_relations, compute_cg_mat
// Calling these functions in the right order is handled by the parent class
// Computations are passed off to the classes FactorBase and RelationGenerator

#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/QuadRelationGenerator.hpp"
#include "ANTL/IndexCalculus/Relation/QuadRelation.hpp"
#include "ANTL/IndexCalculus/FactorBase/QuadFactorBase.hpp"
#include "ANTL/Constants.hpp"

template <class T, class R>
void QuadIndCalc<T,R>::compute_fac_base(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  //cf ANTL/include/ANTL/quadratic/index_calculus/factor_base.hpp
  if ( params.find(Constants::size_fb) == params.end() ) {
    std::cout << "compute_cg_fac_base: value in the first entry is" << params.find(Constants::size_fb)->first << std::endl;
//    NTL::ZZ size_fb(NTL::INIT_VAL, params.find(Constants::size_fb)->first);
//    IndCalc<T,R>::factor_base = QuadFactorBase<T>(size_fb);

    // get factor basis elements using QuadFactorBase
  } else {
    std::cout << "compute_cg_fac_base: params should specify factor basis size as size_fb" << std::endl;
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_relations(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  //cf ANTL/src/quadratic/qo_relation_generator_impl.hpp
//  std::string num_relations_str = "num_relations";
  if ( params.find(Constants::num_relations) == params.end() ) {

    // find relations by calling RelationGenerator's get_relation function
    // Call this function num_relations times

  } else {
    std::cout << "compute_cg_relations: params should specify factor basis size as size_fb" << std::endl;
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_mat(IOrder const &quad_order, std::map<std::string, std::string> const &params) {
  //cf ANTL/include/ANTL/quadratic/index_calculus/qo_relation_matrix.hpp
  IndCalc<T,R>::rels_mat.SetDims(IndCalc<T,R>::relations.size(), IndCalc<T,R>::factor_base.get_size());

  // set matrix using functions in QuadRelation
};

#endif //QUADCLASSGROUPINDCALC_H
