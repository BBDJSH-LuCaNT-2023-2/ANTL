#ifndef RELATION_GENERATOR_H
#define RELATION_GENERATOR_H

#include "ANTL/IndexCalculus/Relation/Relation.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"

namespace ANTL {
  class RelationGenerator {
  public:
    // constructors and destructor
    RelationGenerator();

    virtual ~RelationGenerator() = 0;

    // initialization (set curve and relation generation method)
    virtual long get_relation(Relation &rel);

    RelationGenerator & operator = (const Relation &gen);

  private:
    // Order that relations are generated for
    IOrder *QO = NULL;

    // factor base associated with this relation generator
    FactorBase *FB = NULL;

    std::vector <Relation> rels;
    long total_rels_found;
  };
}

#endif

