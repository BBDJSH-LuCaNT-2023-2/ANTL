#ifndef QUAD_RELATION_GENERATOR_H
#define QUAD_RELATION_GENERATOR_H

#include "ANTL/IndexCalculus/RelationGenerator/RelationGenerator.hpp"
#include "ANTL/IndexCalculus/Relation/QuadRelation.hpp"

namespace ANTL
{

// forward declarations
template <class T> class QuadRelation;

// class
template <class T>
class QuadRelationGenerator : public RelationGenerator
{
  QuadRelationGenerator();
  virtual ~QuadRelationGenerator() = 0;
  // initialization (set curve and relation generation method)
  virtual long get_relation(QuadRelation<T> &rel);
};

} // ANTL

#include "src/IndexCalculus/RelationGenerator/QuadRelationGenerator_impl.hpp"

#endif // guard

