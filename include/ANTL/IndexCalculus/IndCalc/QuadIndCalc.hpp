#ifndef QUADCLASSGROUPINDCALC_H
#define QUADCLASSGROUPINDCALC_H
#include <string>
#include "IndCalc.hpp"
#include "ANTL/IndexCalculus/IndCalc/IndCalc.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/QuadRelationGenerator.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include <NTL/ZZ.h>
#include <vector>

template <class T, class R> // type of Order
class QuadIndCalc:public IndCalc<T,R> {
  /*TODO this class is an example to show which functions a class inheriting from IndCalc needs to implement.
  When we have a fully developed class inheriting from ClassGroupIndCalc this file can be deleted */

public:
  // to make this a concrete class we implement the 4 OrderInvariants
  virtual NTL::ZZ class_number() {return NTL::ZZ(0);};
  virtual std::vector<NTL::ZZ> class_group() {std::vector<NTL::ZZ> cg = {NTL::ZZ(4), NTL::ZZ(5)}; return cg;};
  virtual std::vector<T> unit_group() {std::vector<T> ug = {T()}; return ug;};
  virtual R regulator() {R reg = {R()}; return reg;};

  virtual void compute_fac_base(IOrder const &order, std::map<std::string, std::string> const &params);
  virtual void compute_relations(IOrder const &order, std::map<std::string, std::string> const &params);
  virtual void compute_mat(IOrder const &order, std::map<std::string, std::string> const &params);

  QuadIndCalc<T,R>(IOrder const &quad_order, std::map<std::string, std::string> const &params) : IndCalc<T,R>(quad_order, params) {};
};

#include "src/IndexCalculus/IndCalc/QuadIndCalc_impl.hpp"

#endif //QUADCLASSGROUPINDCALC_H
