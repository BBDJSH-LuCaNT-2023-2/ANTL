#ifndef QUADFACTORBASE_H
#define QUADFACTORBASE_H

#include "ANTL/IndexCalculus/FactorBase/FactorBase.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"

namespace ANTL {

  class QuadFactorBase : public FactorBase {
  public:
    ~QuadFactorBase() {
      std::cout << "desc for QuadFB " << this << std::endl;
    }
//    QuadFactorBase & operator = (const QuadFactorBase &fb);
    using FactorBase::FactorBase; // inherit the construtors
  };
}

#include "src/IndexCalculus/FactorBase/QuadFactorBase_impl.hpp"

#endif //QUADFACTORBASE_H

