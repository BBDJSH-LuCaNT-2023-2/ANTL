#ifndef REGULATOR_LENSTRA_DATA_H
#define REGULATOR_LENSTRA_DATA_H

#include <ANTL/L_function/L_function.hpp>
#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <ANTL/common.hpp>

NTL_CLIENT
using namespace ANTL;

namespace ANTL {
template <class T> class RegulatorLenstraData {
private:
  QuadraticOrder<T> quadratic_order;
  L_function<T> l_function;

public:
  QuadraticOrder<T> get_quadratic_order();
  L_function<T> get_l_function();
};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class T>
QuadraticOrder<T> RegulatorLenstraData<T>::get_quadratic_order() {
  return quadratic_order;
}

template <class T> L_function<T> RegulatorLenstraData<T>::get_l_function() {
  return l_function;
}
} // namespace ANTL
#endif
