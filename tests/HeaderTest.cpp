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
#include <NTL/ZZ.h>
#include "ANTL/IndexCalculus/IndCalc/QuadIndCalc.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Constants.hpp"

using namespace Constants;
using namespace NTL;
using namespace ANTL;

std::map<std::string, std::string> get_params(std::string max_num_tests_str) {
  std::map<std::string, std::string> params {{num_relations, "2"}, {size_fb, "3"}, {bound_fb, "5"},
                                             {max_num_tests, max_num_tests_str}};
  return params;
}

using namespace std;

int main() {

  {
    IOrder order = IOrder();
    long expected_num_relations = 2;
    long expected_size_fb = 3;
    long expected_bound_fb = 5;
    long expected_max_num_tests = 7;
    std::map <std::string, std::string> x;
    x = get_params("7");
    QuadIndCalc <double, double> ind_calc1 = QuadIndCalc<double, double>::create(order, x);

    cout << ind_calc1.get_relation_generator()->get_size_fb() << endl;
    cout << ind_calc1.get_relation_generator()->get_max_num_tests() << endl;
    cout << ind_calc1.get_factor_base()->get_size_fb() << endl;
    cout << ind_calc1.get_factor_base()->get_bound() << endl;
  }
}
