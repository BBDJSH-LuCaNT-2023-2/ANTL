#ifndef RELATION_GENERATOR_H
#define RELATION_GENERATOR_H

#include <string>
#include "ANTL/IndexCalculus/Relation/Relation.hpp"
#include "ANTL/Interface/OrderInvariants.hpp"
#include "ANTL/Constants.hpp"
#include <ANTL/common.hpp>

using namespace ANTL;

namespace ANTL {
  class RelationGenerator {
  public:
    // constructors and destructor
    RelationGenerator(IOrder const &order, std::map<std::string, std::string> const &params, FactorBase const &fb) :
    QO(order), FB(fb) {
      if ( params.find(Constants::size_fb) == params.end() ) {
        std::cout << "RelationGenerator: size_fb should be set" << std::endl;
      } else {
        size_fb = std::stoi(params.find(Constants::size_fb)->second);
      }
    };

    // initialization (set curve and relation generation method)
    virtual long get_relation(Relation &rel) {return 0;}

    RelationGenerator & operator = (const RelationGenerator &gen);

    long get_size_fb() {return size_fb;}

  protected:
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

