#ifndef QUADFACTORBASE_H
#define QUADFACTORBASE_H

#include "ANTL/IndexCalculus/FactorBase/FactorBase.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"

template < class T >
class QuadraticOrder; // an order that inherits from IOrder

namespace ANTL {

  class QuadFactorBase : public FactorBase {
  public:
    using FactorBase::FactorBase; // inherit the construtors
    QuadFactorBase & operator = (const QuadFactorBase &fb);
  };
}

#include "src/IndexCalculus/FactorBase/QuadFactorBase_impl.hpp"

#endif //QUADFACTORBASE_H

