#ifndef QUADCLASSGROUPINDCALC_H
#define QUADCLASSGROUPINDCALC_H
#include <string>
#include "ANTL/IndexCalculus/IndCalc/IndCalc.hpp"

template <class T, class R> // type of Order
class QuadIndCalc:IndCalc<R> {
  /*TODO this class is an example to show which functions a class inheriting from IndCalc needs to implement.
  When we have a fully developed class inheriting from ClassGroupIndCalc this file can be deleted */

public:
  //TODO when QuadraticOrder<T> is done we should replace IOrder with it in this class
  QuadIndCalc<T,R>(IOrder const &quad_order, std::map<std::string, std::string> const &params) : IndCalc<R>(quad_order, params) {IndCalc<R>::setup_cg_mat(quad_order, params);}
  virtual NTL::ZZ class_number() {return NTL::ZZ(0);};
  virtual std::vector<NTL::ZZ> class_group() {std::vector<NTL::ZZ> cg = {NTL::ZZ(4), NTL::ZZ(5)}; return cg;};

  virtual void compute_cg_fac_base(IOrder const &order, std::map<std::string, std::string> const &params);
  virtual void compute_cg_relations(IOrder const &order, std::map<std::string, std::string> const &params);
  virtual void compute_cg_mat(IOrder const &order, std::map<std::string, std::string> const &params);
};

#include "src/IndexCalculus/IndCalc/QuadIndCalc_impl.hpp"

#endif //QUADCLASSGROUPINDCALC_H
