#ifndef FACTORBASE_H
#define FACTORBASE_H

#include <vector>
#include "ANTL/Interface/Multiplicative.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include <ANTL/common.hpp>

using namespace ANTL;

namespace ANTL
{
  class FactorBase {
  protected:
    // Order associated with this factor base
    IOrder *order;// = NULL;

    // size of factor base
    long size_fb; // = 0;

    // bound on factor base primes (degree or integer, depending on base type)
    long bound; // = 0;

    std::vector<IMultiplicative> factor_base;
  public:
    FactorBase(IOrder &new_order) {
      order = &new_order;
    }
    ~FactorBase();

    // accessors
    std::vector<IMultiplicative>& get_fb() {
      return factor_base;
    }

    int get_size() const {
      return size_fb;
    }

    int get_bound() {
      return bound;
    }

    void initialize_size(long size_in);
    void initialize_bound(long bound_in);
  };
}

#endif //FACTORBASE_H
