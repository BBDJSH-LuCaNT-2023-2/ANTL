#ifndef ANTL_VORONOI_METHODS_H
#define ANTL_VORONOI_METHODS_H

#include "CubicIdeal.hpp"


template<typename Type, typename PType>
class CubicIdeal;

template <typename Type, typename PType>
class VoronoiMethods{

public:

/**
* @brief turns the ideal into prepared form.
* @params It is assumed that plattice is the puncture of the coeff_matrix of ideal1
* since this is only ever called by CubicIdeal.make_prepared, this is not an issue.
*/
void make_prepared(CubicIdeal<Type, PType> & ideal1, PType plattice[2][2]);

/**
* @brief virtual stub for converting an ideal's basis to Voronoi form.
* prints an error as the virtual method should not be used.
*/
virtual void make_voronoi_basis(CubicIdeal<Type, PType> & ideal1){
  std::cout << "using the base class method: Error" << std::endl;
};


protected:

// data_members for converting an ideal to a prepared basis
static Type dummy1,
  pm, pm1, pminus, pNext, pBefore,
  qm, qm1, qminus, qNext, qBefore,
  a1,a2,
  m,
  rR;
//            a = reducedIndexForm[0],
//            b = reducedIndexForm[1],
//            c = reducedIndexForm[2],
//            d = reducedIndexForm[3];

static Type temp_lb[3][3];
static PType alpha0, alpha1, alpha2, E;
static PType temp_pb[2][2];

// data members for computing a Voronoi basis
bool omegaDecision[5];
Type omegaMatrix[3][5];

private:




}; // close class def

//define statics
template<typename T, typename PT> T VoronoiMethods<T,PT>::dummy1;
template<typename T, typename PT> T VoronoiMethods<T,PT>::pm;
template<typename T, typename PT> T VoronoiMethods<T,PT>::pm1;
template<typename T, typename PT> T VoronoiMethods<T,PT>::pminus;
template<typename T, typename PT> T VoronoiMethods<T,PT>::pNext;
template<typename T, typename PT> T VoronoiMethods<T,PT>::pBefore;
template<typename T, typename PT> T VoronoiMethods<T,PT>::qm;
template<typename T, typename PT> T VoronoiMethods<T,PT>::qm1;
template<typename T, typename PT> T VoronoiMethods<T,PT>::qminus;
template<typename T, typename PT> T VoronoiMethods<T,PT>::qNext;
template<typename T, typename PT> T VoronoiMethods<T,PT>::qBefore;
template<typename T, typename PT> T VoronoiMethods<T,PT>::a1;
template<typename T, typename PT> T VoronoiMethods<T,PT>::a2;
template<typename T, typename PT> T VoronoiMethods<T,PT>::m;
template<typename T, typename PT> T VoronoiMethods<T,PT>::rR;
//            a = reducedIndexForm[0],
//            b = reducedIndexForm[1],
//            c = reducedIndexForm[2],
//            d = reducedIndexForm[3];

template<typename T, typename PT> T VoronoiMethods<T,PT>::temp_lb[3][3];
template<typename T, typename PT> PT VoronoiMethods<T,PT>::alpha0;
template<typename T, typename PT> PT VoronoiMethods<T,PT>::alpha1;
template<typename T, typename PT> PT VoronoiMethods<T,PT>::alpha2;
template<typename T, typename PT> PT VoronoiMethods<T,PT>::E;
template<typename T, typename PT> PT VoronoiMethods<T,PT>::temp_pb[2][2];

#include "../../../src/Cubic/VoronoiMethods.cpp"



#endif // guard
