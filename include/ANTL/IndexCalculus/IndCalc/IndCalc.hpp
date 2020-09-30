#ifndef CLASSGROUPINDCALC_H
#define CLASSGROUPINDCALC_H

#include <NTL/ZZX.h>
#include <vector>
#include <string>
#include <map>
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Interface/Multiplicative.hpp"
#include "ANTL/IndexCalculus/FactorBase/FactorBase.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/RelationGenerator.hpp"
#include <ANTL/common.hpp>

using namespace ANTL;

template <class R> // type of regulator
class IndCalc {
public:
  IndCalc<R>(IOrder const &order, std::map<std::string, std::string> const &params) {setup_cg_mat(order, params);}

  // subclasses should implement class_number to return the class_number
  virtual NTL::ZZ class_number();
  // subclasses should implement class_group to return the class_group
  virtual std::vector<NTL::ZZ> class_group();

  FactorBase cg_factor_base;
  RelationGenerator cg_relations_generator;
  NTL::Mat<long> cg_mat;

  void setup_cg_fac_base(IOrder const &order, std::map<std::string, std::string> const &params);
  void setup_cg_relations(IOrder const &order, std::map<std::string, std::string> const &params);
  void setup_cg_mat(IOrder const &order, std::map<std::string, std::string> const &params);

protected:
  // subclasses should implement compute_cg_fac_base to set the cg_factor_base variable
  virtual void compute_cg_fac_base(IOrder const &order, std::map<std::string, std::string> const &params) {};
  // subclasses should implement compute_cg_relations to set the cg_relations variable
  virtual void compute_cg_relations(IOrder const &order, std::map<std::string, std::string> const &params) {};
  // subclasses should implement compute_cg_mat to set the cg_mat variable
  virtual void compute_cg_mat(IOrder const &order, std::map<std::string, std::string> const &params) {};
};

template <class R> void IndCalc<R>::setup_cg_fac_base(const IOrder &order, const std::map<std::string, std::string> &params) {
  compute_cg_fac_base(order, params);
  //TODO
//  if (cg_factor_base.size() == 0) {
    // raise error
//  }
}

template <class R> void IndCalc<R>::setup_cg_relations(IOrder const &order, std::map<std::string, std::string> const &params) {
  setup_cg_fac_base(order, params);
  compute_cg_relations(order, params);
  //TODO
//  if (cg_relations.size() == 0) {
    // raise error
//  }
}

template <class R> void IndCalc<R>::setup_cg_mat(IOrder const &order, std::map<std::string, std::string> const &params) {
  setup_cg_relations(order, params);
  compute_cg_mat(order, params);
  if (cg_mat.NumRows() == 0 or cg_mat.NumCols() == 0) {
    // raise error
  }
}
#endif //CLASSGROUPINDCALC_H
