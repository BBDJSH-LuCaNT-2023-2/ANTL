#ifndef QUADFACTORBASE_H
#define QUADFACTORBASE_H

#include <ANTL/IndexCalculus/FactorBase/FactorBase.hpp>
#include <ANTL/Interface/OrderInvariants.hpp>

namespace ANTL {

  class QuadFactorBase : public FactorBase {
  public:
    using FactorBase::FactorBase; // inherit the construtors
  };
}

#include "src/IndexCalculus/FactorBase/QuadFactorBase_impl.hpp"

#endif //QUADFACTORBASE_H

