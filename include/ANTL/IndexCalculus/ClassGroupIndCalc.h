//
// Created by David Marquis on 2020-09-16.
//

#ifndef CLASSGROUPINDCALC_H
#define CLASSGROUPINDCALC_H

#include <NTL/ZZX.h>
#include <vector>
#include <string>
#include <map>
#include "../Interface/OrderInvariants.h"
#include "../Interface/Multiplicative.h"

class ClassGroupIndCalc {
public:
  ClassGroupIndCalc(IOrder const &order, std::map<std::string, std::string> const &params) {setup_cg_mat(order, params);}

  // subclasses should implement class_number to return the class_number using cg_mat
  virtual NTL::ZZ class_number() {};
  // subclasses should implement class_group to return the class_group using cg_mat
  virtual std::string class_group() {};

  std::vector<IMultiplicative> cg_factor_base;
  std::vector<IMultiplicative> cg_relations;
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

void ClassGroupIndCalc::setup_cg_fac_base(const IOrder &order, const std::map<std::string, std::string> &params) {
  compute_cg_fac_base(order, params);
  if (cg_factor_base.size() == 0) {
    // raise error
  }
}

void ClassGroupIndCalc::setup_cg_relations(IOrder const &order, std::map<std::string, std::string> const &params) {
  setup_cg_fac_base(order, params);
  compute_cg_relations(order, params);
  if (cg_relations.size() == 0) {
    // raise error
  }
}

void ClassGroupIndCalc::setup_cg_mat(IOrder const &order, std::map<std::string, std::string> const &params) {
  setup_cg_relations(order, params);
  compute_cg_mat(order, params);
  if (cg_mat.NumRows() == 0 or cg_mat.NumCols() == 0) {
    // raise error
  }
}
#endif //CLASSGROUPINDCALC_H
