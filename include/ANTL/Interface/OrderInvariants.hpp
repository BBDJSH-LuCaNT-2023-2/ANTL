//
// Created by David Marquis on 2020-09-17.
//

#ifndef ORDERINVARIANTS_H
#define ORDERINVARIANTS_H

#include <NTL/ZZ.h>
#include <vector>

class IOrder {
public:
  // interface for objects that behave like orders in algebraic number fields
  // and have invariants like class number and regulator

  IOrder &operator = (IOrder const &order);
  ~IOrder() {
    std::cout << "desc for IOrder " << this << std::endl;
  }

};

class IClassGroup:public IOrder {
public:
  virtual std::vector<NTL::ZZ> class_group() = 0;
};

class IClassNumber:public IOrder {
public:
  virtual NTL::ZZ class_number() = 0;
};

template <class T> // type of a unit
class IUnitGroup:public IOrder {
public:
  virtual std::vector<T> unit_group() = 0;
};

template <class R> // type of the regulator
class IRegulator:public IOrder {
public:
  virtual R regulator() = 0;
};

#endif //ORDERINVARIANTS_H
