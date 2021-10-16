// This file is used to test that header files can be built
// Currently this only tests a subset of the headers. We can grow this set as we go

#include "ANTL/IndexCalculus/IndCalc/IndCalc.hpp"
#include "ANTL/IndexCalculus/IndCalc/QuadIndCalc.hpp"

#include "ANTL/Interface/Multiplicative.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/IndexCalculus/Relation/Relation.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/RelationGenerator.hpp"

#include "ANTL/IndexCalculus/FactorBase/QuadFactorBase.hpp"
#include "ANTL/IndexCalculus/Relation/QuadRelation.hpp"
#include "ANTL/IndexCalculus/RelationGenerator/QuadRelationGenerator.hpp"

#include <string>
#include <map>
#include <iostream>
#include <typeinfo>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include "ANTL/IndexCalculus/IndCalc/QuadIndCalc.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Quadratic/QuadraticOrder.hpp"
#include "ANTL/Constants.hpp"

using namespace Constants;
using namespace NTL;
using namespace ANTL;
using namespace std;

std::map<std::string, std::string> get_params(std::string max_num_tests_str) {
  std::map<std::string, std::string> params {{num_relations, "2"}, {size_fb, "3"}, {bound_fb, "5"},
                                             {max_num_tests, max_num_tests_str}};
  return params;
}

int main() {

  QuadraticOrder<ZZ> order = QuadraticOrder<ZZ>(ZZ(13));
  {

//    IOrder order = IOrder();
    long expected_num_relations = 2;
    long expected_size_fb = 3;
    long expected_bound_fb = 5;
    long expected_max_num_tests = 7;
    std::map <std::string, std::string> x;
    x = get_params("7");
    auto ind_calc = QuadIndCalc<ZZ, RR>::create(order, x);

//    auto ind_calc = new QuadIndCalc<double,double>();
//    ind_calc->relation_generator = std::move(reln_generator);
//    ind_calc->factor_base = std::move(fac_base);
//    ind_calc->setup_mat();

//    cout << ind_calc->get_relation_generator()->get_size_fb() << endl;
//    cout << ind_calc->get_relation_generator()->get_max_num_tests() << endl;
//    cout << ind_calc->get_factor_base()->get_size_fb() << endl;
//    cout << ind_calc->get_factor_base()->get_bound() << endl;
//    cout << "I should be able to call the destructor on ind_calc" << endl;

//    ind_calc.reset(); // this by itself is ok
//    ind_calc->get_relation_generator().reset(); // this by itself is ok
//    ind_calc->get_relation_generator().reset(); // this by itself is ok
  }
  cout << "ind calc destroyed but order still alive" << endl;
//  cout << &order << endl;
}
