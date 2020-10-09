#ifndef QUAD_RELATION_GENERATOR_H
#define QUAD_RELATION_GENERATOR_H

#include "ANTL/IndexCalculus/RelationGenerator/RelationGenerator.hpp"
#include "ANTL/IndexCalculus/Relation/QuadRelation.hpp"
#include "ANTL/Constants.hpp"

namespace ANTL
{

// forward declarations
template <class T> class QuadRelation;

// class
template <class T>
class QuadRelationGenerator : public RelationGenerator
{
public:
  QuadRelationGenerator<T>(IOrder const &order, std::map<std::string, std::string> const &params, FactorBase const &fb) :
  RelationGenerator(order, params, fb) {};
  QuadRelationGenerator & operator = (const QuadRelationGenerator &fb);
};

} // ANTL

#include "src/IndexCalculus/RelationGenerator/QuadRelationGenerator_impl.hpp"

#endif // guard

