#ifndef ANTL_BASIC_VORONOI_HPP
#define ANTL_BASIC_VORONOI_HPP

#include "FundUnitStrategy.hpp"


template<typename Type, typename PType>
class CubicOrder;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class FundUnitStrategy;
template<typename Type, typename PType>
class BasicVoronoi : public FundUnitStrategy<Type, PType> {

public:

  BasicVoronoi();

  ~BasicVoronoi(){
    //delete mul_holder;

  };

void compute(std::vector<CubicElement<Type, PType>> & unitvec, CubicOrder<Type, PType> * ord, bool realorder){
  if (realorder == true){
    fundamental_unit_real(unitvec, ord);
  }
  else{
    fundamental_unit_complex(unitvec, ord);
  }
}


protected:
  void fundamental_unit_complex(std::vector<CubicElement<Type, PType>> & unitvec, CubicOrder<Type, PType> * ord);

  void fundamental_unit_real(std::vector<CubicElement<Type, PType>> & unitvec, CubicOrder<Type, PType> * ord);


private:

  bool completeCycle;

  CubicIdeal<Type, PType> L1 = CubicIdeal<Type, PType>(NULL);
  CubicIdeal<Type, PType> L2 = CubicIdeal<Type, PType>(NULL);;

  CubicElement<Type, PType> epsilon1 = CubicElement<Type, PType>(NULL);
  CubicElement<Type, PType> adj_minimum = CubicElement<Type, PType>(NULL);
  CubicElement<Type, PType> epsilon2 = CubicElement<Type, PType>(NULL);


// vectors used in the real fundamental unit calculations
std::vector<CubicIdeal<Type, PType>> x_cycle;         // container to hold x-cycle
std::vector<CubicElement<Type, PType>> adj_minima_vec;    // container for holding adjacent min



}; //end class definition
#include "../../../../src/Cubic/FundamentalUnits/BasicVoronoi.cpp"



#endif
