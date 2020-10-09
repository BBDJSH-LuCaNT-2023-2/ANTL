#ifndef FACTORBASE_H
#define FACTORBASE_H

#include <vector>
#include <iostream>
#include <string>
#include "ANTL/Interface/Multiplicative.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Constants.hpp"
#include <ANTL/common.hpp>

using namespace ANTL;

namespace ANTL
{
  class FactorBase {
  protected:
    // Order associated with this factor base
    IOrder const &order;

    // size of factor base
    long size_fb;

    // bound on factor base primes (degree or integer, depending on base type)
    long bound;

    std::vector<IMultiplicative> factor_base;
  public:
    long get_size_fb() {return size_fb;}
    FactorBase & operator = (const FactorBase &fb);

    FactorBase(IOrder const &new_order, std::map<std::string, std::string> const &params) : order(new_order) {
      if ( params.find(Constants::size_fb) == params.end() ) {
        std::cout << "FactorBase: size_fb should be set" << std::endl;
      } else {
        size_fb = std::stoi(params.find(Constants::size_fb)->second);
      }
      if ( params.find(Constants::bound_fb) == params.end() ) {
        std::cout << "FactorBase: bound_fb should be set" << std::endl;
      } else {
        bound = std::stoi(params.find(Constants::bound_fb)->second);
      }
    }

    // accessors
    std::vector<IMultiplicative>& get_fb() {
      return factor_base;
    }

    long get_size() const {
      return size_fb;
    }

    long get_bound() {
      return bound;
    }

    void initialize_size(long size_in);
    void initialize_bound(long bound_in);
  };
}

#endif //FACTORBASE_H
