//
// Created by David Marquis on 2020-09-16.
//

#ifndef RUNNTL_QUADCLASSGROUPINDCALC_H
#define RUNNTL_QUADCLASSGROUPINDCALC_H
#include <string>
#include "ClassGroupIndCalc.h"
#include "Order.h"

class QuadClassGroupIndCalc:ClassGroupIndCalc {
public:
  QuadClassGroupIndCalc(IOrder const &nf_order, std::map<std::string, std::string> const &params) : ClassGroupIndCalc(nf_order, params) {ClassGroupIndCalc::setup_cg_mat(nf_order, params);}
  virtual void compute_cg_fac_base(IOrder const &quad_order, std::map<std::string, std::string> const &params) {};
  virtual void compute_cg_relations(IOrder const &quad_order, std::map<std::string, std::string> const &params) {};
  virtual void compute_cg_mat(IOrder const &quad_order, std::map<std::string, std::string> const &params) {};
  virtual NTL::ZZ class_number() {return NTL::ZZ(0);};
  virtual std::string class_group() {std::string cg = "unimplemented"; return cg;};

};

#endif //RUNNTL_QUADCLASSGROUPINDCALC_H
