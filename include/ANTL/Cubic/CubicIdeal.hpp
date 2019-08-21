#ifndef ANTL_CUBIC_IDEAL_H
#define ANTL_CUBIC_IDEAL_H

template<typename Type,typename PType>
class CubicIdeal{

public:


protected:

// holds a reference to the order in which the ideal sits in
CubicOrder<Type, PType> * my_order;

// a Z basis for the ideal, should be of size 3
// Strongly consider creating a Basis class
Type * my_basis;
Type denom;


private:


}

#endif
