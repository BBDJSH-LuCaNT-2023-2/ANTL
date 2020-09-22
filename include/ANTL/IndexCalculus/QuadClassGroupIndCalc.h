//
// Created by David Marquis on 2020-09-16.
//

#ifndef QUADCLASSGROUPINDCALC_H
#define QUADCLASSGROUPINDCALC_H
#include <string>
#include "ClassGroupIndCalc.h"

class QuadClassGroupIndCalc:ClassGroupIndCalc {
  /*TODO this class is an example to show which functions a class inheriting from ClassGroupIndCalc needs to implement.
  When we have a fully developed class inheriting from ClassGroupIndCalc this file can be deleted */
public:
  QuadClassGroupIndCalc(IOrder const &quad_order, std::map<std::string, std::string> const &params) : ClassGroupIndCalc(quad_order, params) {ClassGroupIndCalc::setup_cg_mat(quad_order, params);}
  virtual void compute_cg_fac_base(IOrder const &quad_order, std::map<std::string, std::string> const &params) {};
  virtual void compute_cg_relations(IOrder const &quad_order, std::map<std::string, std::string> const &params) {};
  virtual void compute_cg_mat(IOrder const &quad_order, std::map<std::string, std::string> const &params) {};
  virtual NTL::ZZ class_number() {return NTL::ZZ(0);};
  virtual std::string class_group() {std::string cg = "unimplemented"; return cg;};
};

#endif //QUADCLASSGROUPINDCALC_H
