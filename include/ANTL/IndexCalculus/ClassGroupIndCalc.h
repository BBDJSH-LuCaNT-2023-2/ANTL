//
// Created by David Marquis on 2020-09-16.
//

#ifndef RUNNTL_CLASSGROUPINDCALC_H
#define RUNNTL_CLASSGROUPINDCALC_H

#include "OrderInvariants.h"
#include <NTL/ZZX.h>
#include <vector>

#include <string>
#include <map>
#include "Elem.h"

class ClassGroupIndCalc {
public:
  ClassGroupIndCalc(IOrder const &nf_order, std::map<std::string, std::string> const &params) {setup_cg_mat(nf_order, params);}

  virtual NTL::ZZ class_number() {};
  virtual std::string class_group() {};

  std::vector<IMultiplicative> cg_relations;
  std::vector<IMultiplicative> cg_factor_base;
  NTL::Mat<long> cg_mat;

  void setup_cg_fac_base(IOrder const &nf_order, std::map<std::string, std::string> const &params);
  void setup_cg_relations(IOrder const &nf_order, std::map<std::string, std::string> const &params);
  void setup_cg_mat(IOrder const &nf_order, std::map<std::string, std::string> const &params);

protected:
  virtual void compute_cg_fac_base(IOrder const &nf_order, std::map<std::string, std::string> const &params) {};
  virtual void compute_cg_relations(IOrder const &nf_order, std::map<std::string, std::string> const &params) {};
  virtual void compute_cg_mat(IOrder const &nf_order, std::map<std::string, std::string> const &params) {};
};

void ClassGroupIndCalc::setup_cg_fac_base(const IOrder &nf_order, const std::map<std::string, std::string> &params) {
  compute_cg_fac_base(nf_order, params);
  if (cg_factor_base.size() == 0) {
    // raise error
  }
}

void ClassGroupIndCalc::setup_cg_relations(IOrder const &nf_order, std::map<std::string, std::string> const &params) {
  setup_cg_fac_base(nf_order, params);
  compute_cg_relations(nf_order, params);
  if (cg_relations.size() == 0) {
    // raise error
  }
}

void ClassGroupIndCalc::setup_cg_mat(IOrder const &nf_order, std::map<std::string, std::string> const &params) {
  setup_cg_relations(nf_order, params);
  compute_cg_mat(nf_order, params);
  if (cg_mat.NumRows() == 0 or cg_mat.NumCols() == 0) {
    // raise error
  }
}
#endif //RUNNTL_CLASSGROUPINDCALC_H
