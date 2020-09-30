#ifndef RELATION_GENERATOR_H
#define RELATION_GENERATOR_H

#include "ANTL/IndexCalculus/Relation/Relation.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include <ANTL/common.hpp>

using namespace ANTL;

namespace ANTL {
  class RelationGenerator {
  public:
    // constructors and destructor
    RelationGenerator();
    RelationGenerator(IOrder const &order, std::map<std::string, std::string> const &params, FactorBase const &fb) :
    QO(order), FB(fb) {size_fb = std::stoi(params.find("size_fb")->first);};

    virtual ~RelationGenerator() = 0;

    // initialization (set curve and relation generation method)
    virtual long get_relation(Relation &rel);

    RelationGenerator & operator = (const Relation &gen);

  private:
    // size of factor base
    long size_fb;

    // Order that relations are generated for
    IOrder QO;

    // factor base associated with this relation generator
    FactorBase FB;

    std::vector <Relation> rels;
    long total_rels_found;
  };
}

#endif

