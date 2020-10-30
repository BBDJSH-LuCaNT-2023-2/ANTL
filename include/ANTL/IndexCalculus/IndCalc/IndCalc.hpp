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
class IndCalc: IClassGroup, IClassNumber, IUnitGroup<T>, IRegulator<R> {
public:
  std::vector<Relation> relations;
  NTL::Mat<ZZ> rels_mat;

  // subclasses must implement class_number, class group, unit group, regulator
  NTL::ZZ class_number() {};
  std::vector<NTL::ZZ> class_group() {}
  std::vector<T> unit_group() {}
  R regulator() {}

  virtual FactorBase* get_factor_base() {return factor_base.get();};
  virtual RelationGenerator* get_relation_generator() {return relation_generator.get();};

protected:
  void setup_fac_base(IOrder const &order, std::map<std::string, std::string> const &params);
  void setup_relations(IOrder const &order, std::map<std::string, std::string> const &params);
  void setup_mat(IOrder const &order, std::map<std::string, std::string> const &params);

  // subclasses should implement compute_fac_base to set the factor_base variable
  virtual void compute_fac_base(IOrder const &order, std::map<std::string, std::string> const &params) {};
  // subclasses should implement compute_relations to set the cg_relations variable
  virtual void compute_relations(IOrder const &order, std::map<std::string, std::string> const &params) {};
  // subclasses should implement compute_mat to set the cg_mat variable
  virtual void compute_mat(IOrder const &order, std::map<std::string, std::string> const &params) {};

  IndCalc<T,R>() = default;
  IndCalc<T,R>(std::shared_ptr<FactorBase> factor_base, std::shared_ptr<RelationGenerator> relation_generator) :
    factor_base(factor_base),
    relation_generator(relation_generator) {}

private:
  std::shared_ptr<FactorBase> factor_base;
  std::shared_ptr<RelationGenerator> relation_generator;
};

template <class T, class R> void IndCalc<T,R>::setup_fac_base(const IOrder &order, const std::map<std::string, std::string> &params) {
  std::cout << "in setup_fac_base" << std::endl;
  compute_fac_base(order, params);
  //TODO
//  if (factor_base.size() == 0) {
    // raise error
//  }
}

template <class T, class R> void IndCalc<T,R>::setup_relations(IOrder const &order, std::map<std::string, std::string> const &params) {
  setup_fac_base(order, params);
  std::cout << "in setup_relations" << std::endl;
  compute_relations(order, params);
  //TODO
//  if (relations.size() == 0) {
    // raise error
//  }
}

template <class T, class R> void IndCalc<T,R>::setup_mat(IOrder const &order, std::map<std::string, std::string> const &params) {
  setup_relations(order, params);
  compute_mat(order, params);
  if (rels_mat.NumRows() == 0 or rels_mat.NumCols() == 0) {
    // raise error
  }
}
#endif //CLASSGROUPINDCALC_H
