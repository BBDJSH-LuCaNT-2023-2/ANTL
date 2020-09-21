//
// Created by David Marquis on 2020-09-17.
//

#ifndef RUNNTL_ORDERINVARIANTS_H
#define RUNNTL_ORDERINVARIANTS_H

#include <NTL/ZZ.h>

class IOrder {
  // interface for objects that behave like orders in algebraic number fields
  // and have invariants like class_number and regulator
};

class IClassGroup:public IOrder {
public:
  virtual std::string class_group() = 0;
};

class IClassNumber:public IOrder {
public:
  virtual NTL::ZZ class_number() = 0;
};

template <class T> // type of a unit
class IUnitGroup:public IOrder {
public:
  virtual T unit_group() = 0;
};

template <class R> // class for the type of real numbers
class IRegulator:public IOrder {
public:
  virtual NTL::ZZ regulator() = 0;
};

#endif //RUNNTL_ORDERINVARIANTS_H
