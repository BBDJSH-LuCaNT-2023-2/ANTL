#ifndef ANTL_VORONOI_REAL_H
#define ANTL_VORONOI_REAL_H

#include "VoronoiMethods.hpp"


template<typename Type, typename PType>
class VoronoiMethods;


template <typename Type, typename PType>
class VoronoiReal : public VoronoiMethods<Type, PType>{


public:

void make_voronoi_basis(CubicIdeal<Type, PType> & ideal1, bool reduced = true);
protected:

int p1, p2;
CubicOrderReal<Type, PType>  * real_ord;
private:





};




#include "../../../src/Cubic/VoronoiReal.cpp"

#endif // guard
