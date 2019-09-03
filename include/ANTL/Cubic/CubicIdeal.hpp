#ifndef ANTL_CUBIC_IDEAL_H
#define ANTL_CUBIC_IDEAL_H

#include "CubicOrderNF.hpp"
#include "CubicElementNF.hpp"


template<typename Type,typename PType>
class CubicIdeal{

public:


/************** Accessors **********************/
inline CubicElement<Type, PType> * get_gen1() const {return gen1;}
inline CubicElement<Type, PType> * get_gen2() const {return gen2;}
inline CubicElement<Type, PType> * get_gen3() const {return gen3;}
inline Type get_denom () const {return denom;}
inline const CubicOrder<Type, PType> * get_order() const {return my_order;}

Type get_norm();

bool is_principal();
bool is_prime();
bool is_equivalent(const CubicIdeal<Type, PType> & B);
bool is_canonical();



void reduce();
void become_canonical();
void become_prepared();
void voronoi();

/******************** friends *********************/
friend void mul<Type, PType > (CubicIdeal<Type, PType> & A, const CubicIdeal<Type, PType> & B, const CubicIdeal<Type, PType> & C);
protected:


// holds a reference to the order in which the ideal sits in
CubicOrder<Type, PType> * my_order;

// a Z basis for the ideal, should be of size 3
// Strongly consider creating a Basis class
CubicElement<Type, PType> gen1;
CubicElement<Type, PType> gen2;
CubicElement<Type, PType> gen3;
Type denom;

private:

void normalize() = 0;


}; // close class definition

#endif
