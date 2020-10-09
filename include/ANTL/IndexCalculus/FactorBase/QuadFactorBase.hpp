#ifndef QUADFACTORBASE_H
#define QUADFACTORBASE_H

#include "ANTL/IndexCalculus/FactorBase/FactorBase.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"

namespace ANTL {
  // this fake QuadraticOrder is used as a placeholder while ANTL/Quadratic/QuadraticOrder is under development
  template < class T >
  class QuadraticOrder: public IOrder {};

  template <class T> // type of ideals in the order
  class QuadFactorBase : public FactorBase {
  public:

    QuadFactorBase(QuadraticOrder<T> const &new_order, std::map<std::string, std::string> const &params) : FactorBase(new_order, params) {

    }
    QuadFactorBase & operator = (const QuadFactorBase &fb);

  };
}

#include "src/IndexCalculus/FactorBase/QuadFactorBase_impl.hpp"

#endif //QUADFACTORBASE_H

