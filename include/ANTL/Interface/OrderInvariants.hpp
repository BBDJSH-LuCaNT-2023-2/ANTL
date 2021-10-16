//
// Created by David Marquis on 2020-09-17.
//

#ifndef ORDERINVARIANTS_H
#define ORDERINVARIANTS_H

#include <NTL/ZZ.h>
#include <vector>

//template<class T>
//class QuadraticNumber
template <class T, class R> // type of unit, type of regulator
class IOrder {
public:
  // interface for objects that behave like orders in algebraic number fields
  // and have invariants like class number and regulator

//  IOrder() = default;
//  IOrder(const IOrder& rhs) = default;
//  IOrder &operator = (IOrder const &order);
  ~IOrder() {
    std::cout << "desc for IOrder " << this << std::endl;
  }

  // subclasses may implement class_group, class_number, unit_group, regulator
  virtual std::vector<NTL::ZZ> class_group() {
	std::vector<NTL::ZZ> classGroup = {NTL::ZZ(0)};
	return classGroup;
  }

  virtual NTL::ZZ class_number() {
	NTL::ZZ classNumber = NTL::ZZ(0);
	return classNumber;
  }

  virtual std::vector<T> unit_group() {
	T object1, object2;
	std::vector<T> unitGroup = {object1, object2};
	return unitGroup;
  }

  virtual R regulator() {
	R regulator;
	return regulator;
  }

};

#endif //ORDERINVARIANTS_H
