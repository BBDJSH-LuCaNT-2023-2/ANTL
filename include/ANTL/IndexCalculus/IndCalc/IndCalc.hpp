#ifndef CLASSGROUPINDCALC_H
#define CLASSGROUPINDCALC_H

#include <NTL/ZZX.h>
#include <vector>
#include <string>
#include <map>
#include <memory>
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Interface/Multiplicative.hpp"
#include "ANTL/IndexCalculus/FactorBase/FactorBase.hpp"
#include "ANTL/IndexCalculus/Relation/Relation.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/RelationGenerator.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include <ANTL/common.hpp>

using namespace ANTL;

template <class T, class R> // type of unit, type of regulator
class IndCalc : IOrder<T,R> {
public:
  std::vector<Relation> relations;
  NTL::Mat<ZZ> rels_mat;

  // subclasses must implement class_number, class group, unit group, regulator
  NTL::ZZ class_number() {};
  std::vector<NTL::ZZ> class_group() {}
  std::vector<T> unit_group() {}
  R regulator() {}

  virtual FactorBase* const get_factor_base() = 0;
  virtual RelationGenerator* const get_relation_generator() = 0;

//TODO protected:
  void setup_fac_base();
  void setup_relations();
  void setup_mat();

  // subclasses should implement compute_fac_base
  virtual void compute_fac_base() {};
  // subclasses should implement compute_relations
  virtual void compute_relations() {};
  // subclasses should implement compute_mat
  virtual void compute_mat() {};

  IndCalc<T,R>() = default;
};

template <class T, class R> void IndCalc<T,R>::setup_fac_base() {
  compute_fac_base();
  //TODO
//  if (factor_base.size() == 0) {
    // raise error
//  }
}

template <class T, class R> void IndCalc<T,R>::setup_relations() {
  setup_fac_base();
  compute_relations();
  //TODO
//  if (relations.size() == 0) {
    // raise error
//  }
}

template <class T, class R> void IndCalc<T,R>::setup_mat() {
  setup_relations();
  compute_mat();
  if (rels_mat.NumRows() == 0 or rels_mat.NumCols() == 0) {
    // raise error
  }
}
#endif //CLASSGROUPINDCALC_H
