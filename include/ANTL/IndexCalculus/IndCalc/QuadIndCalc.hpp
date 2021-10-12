#ifndef QUADCLASSGROUPINDCALC_H
#define QUADCLASSGROUPINDCALC_H
#include <string>
#include "IndCalc.hpp"
#include "ANTL/IndexCalculus/IndCalc/IndCalc.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include <NTL/ZZ.h>
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/QuadRelationGenerator.hpp"
#include "ANTL/IndexCalculus/Relation/Relation.hpp"
#include "ANTL/IndexCalculus/Relation/QuadRelation.hpp"
#include "ANTL/IndexCalculus/FactorBase/QuadFactorBase.hpp"
#include "ANTL/Constants.hpp"
#include <vector>
#include <memory>

template <class T, class R> // type of Order
class QuadIndCalc:public IndCalc<T,R> {
  /*TODO this class is an example to show which functions a class inheriting from IndCalc needs to implement.
  When we have a fully developed class inheriting from ClassGroupIndCalc this file can be deleted */

public:
  using IndCalc<T,R>::IndCalc; // inherit the constructors
  friend void FactorBase::push_to_fb(IMultiplicative &fb_elem); // ind calc can add stuff to the factor base

  //TODO fill in these stubs after we have implemented an index calculus algorithm for this class
  virtual NTL::ZZ class_number() {return NTL::ZZ(0);};
  virtual std::vector<NTL::ZZ> class_group() {std::vector<NTL::ZZ> cg = {NTL::ZZ(4), NTL::ZZ(5)}; return cg;};
  virtual std::vector<T> unit_group() {std::vector<T> ug = {T()}; return ug;};
  virtual R regulator() {R reg = {R()}; return reg;};

  QuadIndCalc<T,R>() = default;
  QuadIndCalc(const QuadIndCalc&) = delete;
  QuadIndCalc& operator=(const QuadIndCalc&) = delete;
  ~QuadIndCalc() = default;
  QuadIndCalc<T,R>(std::unique_ptr<QuadFactorBase> factor_base, std::unique_ptr<QuadRelationGenerator> relation_generator) :
    factor_base(std::move(factor_base)),
    relation_generator(std::move(relation_generator)) {}

  static std::unique_ptr<QuadIndCalc<ZZ,RR>> create(IOrder<ZZ,RR> const &order, std::map<std::string, std::string> const &params) {
    // two-phase initialization
    // TODO: make this function general

    unique_ptr<QuadFactorBase> fac_base{ new QuadFactorBase(std::move(order), params) };
    unique_ptr<QuadRelationGenerator> reln_generator{ new QuadRelationGenerator(std::move(order), params, fac_base.get())};
    unique_ptr<QuadIndCalc<ZZ,RR>> ind_calc {new QuadIndCalc<ZZ,RR>(std::move(fac_base), std::move(reln_generator))};
    ind_calc->setup_mat();
    return ind_calc;
  }

  virtual RelationGenerator* get_relation_generator() override {return relation_generator.get();};
  virtual FactorBase* get_factor_base() override {return factor_base.get();};

  void compute_fac_base() override;
  void compute_relations() override;
  void compute_mat() override;

private:
  std::unique_ptr<QuadFactorBase> factor_base;
  std::unique_ptr<QuadRelationGenerator> relation_generator;
};

// For now these 3 functions stubs. The code in them is just to suggest what they need to do

template <class T, class R>
void QuadIndCalc<T,R>::compute_fac_base() {
  IMultiplicative fb_elem = IMultiplicative();

  for(int i=0; i < factor_base->get_size(); i++) {
    factor_base->push_to_fb(fb_elem);
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_relations() {
  bool result;
  /* compute_relations does sieving or random exponents to construct a relation */

  Relation quad_relation = Relation();
  for(long i=0; i < get_relation_generator()->get_max_num_tests(); i++) {
    result = get_relation_generator()->get_relation(quad_relation, i);
    if (result) {
      IndCalc<T,R>::relations.push_back(quad_relation);
    }
  }
};

template <class T, class R>
void QuadIndCalc<T,R>::compute_mat() {
  //cf ANTL/include/ANTL/quadratic/index_calculus/qo_relation_matrix.hpp
  IndCalc<T,R>::rels_mat.SetDims(IndCalc<T,R>::relations.size(), factor_base->get_size());

  // set matrix using functions in QuadRelation
};

#include "src/IndexCalculus/IndCalc/QuadIndCalc_impl.hpp"

#endif //QUADCLASSGROUPINDCALC_H
