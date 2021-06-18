#ifndef ANTL_BSGS_VORONOI_HPP
#define ANTL_BSGS_VORONOI_HPP

#include "FundUnitStrategy.hpp"
#include <iterator>
#include <unordered_map>

template<typename Type, typename PType>
class CubicOrder;
template<typename Type, typename PType>
class Cubic;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class FundUnitStrategy;
template<typename Type, typename PType>
class BasicVoronoi;
template<typename Type, typename PType>
class BSGSVoronoi : public FundUnitStrategy<Type, PType> {

public:

  BSGSVoronoi();

  ~BSGSVoronoi(){
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

  void doubling(CubicIdeal<Type, PType> & R, PType & logval );

  void initialize_giant_baby(CubicOrder<Type, PType> * ord);

  void giant_step(CubicIdeal<Type, PType> & currentIdeal, CubicElement<Type, PType> & currentMinimum, PType & logval );
private:

  bool collision, completeCycle;
  Type norm_holder;
  PType distance, realtemp;
  PType giantBound, giant_logarithm;
  ZZHash hash_function;

  CubicIdeal<Type, PType> L1 = CubicIdeal<Type, PType>(NULL);
  CubicIdeal<Type, PType> L2 = CubicIdeal<Type, PType>(NULL);

  CubicElement<Type, PType> epsilon1 = CubicElement<Type, PType>(NULL);
  CubicElement<Type, PType> adj_minimum = CubicElement<Type, PType>(NULL);
  CubicElement<Type, PType> epsilon2 = CubicElement<Type, PType>(NULL);


// vectors used in the real fundamental unit calculations
std::vector<CubicIdeal<Type, PType>> x_cycle;         // container to hold x-cycle
std::vector<CubicElement<Type, PType>> adj_minima_vec;    // container for holding adjacent min


//static std::unordered_multimap<Type, CubicIdeal<Type,PType> > babysteps;
CubicIdeal<Type, PType> bigFoot = CubicIdeal<Type, PType>(NULL);
CubicElement<Type, PType> giant_min = CubicElement<Type, PType>(NULL);


std::unordered_multimap<Type, std::size_t, ZZHash, ZZEqual > babysteps = std::unordered_multimap<Type, std::size_t, ZZHash, ZZEqual>( 1, hash_function, ZZEqual() );
}; //end class definition



#include "../../../../src/Cubic/FundamentalUnits/BSGSVoronoi.cpp"



#endif
