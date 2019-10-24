#ifndef ANTL_VORONOI_COMPLEX_H
#define ANTL_VORONOI_COMPLEX_H

#include "VoronoiMethods.hpp"


template<typename Type, typename PType>
class VoronoiMethods;


template <typename Type, typename PType>
class VoronoiComplex : public VoronoiMethods<Type, PType>{


public:

void make_voronoi_basis(CubicIdeal<Type, PType> & ideal1, bool reduced=true);
protected:

CubicElement<Type, PType> placeholder;

private:





};




#include "../../../src/Cubic/VoronoiComplex.cpp"

#endif // guard
