#ifndef ANTL_CUBIC_IDEAL_H
#define ANTL_CUBIC_IDEAL_H

/**
 * @file CubicIdeal.hpp
 * @author Randy Yee
 * @remarks Class representing a generic cubic ideal.
 */


#include "CubicOrder.hpp"
#include "CubicElement.hpp"
#include "CubicOrderReal.hpp"
#include "../Arithmetic/QQ.hpp"

using namespace ANTL;

//forward declarations
template<typename Type, typename PType>
class CubicOrder;
template<typename Type,typename PType>
class CubicOrderReal;
template<typename Type, typename PType>
class CubicElement;
template<typename Type, typename PType>
class MultiplyStrategyWilliams;
template<typename Type,typename PType>
class IdealMultiplicationStrategy;
template<typename Type, typename PType>
class VoronoiMethods;

template<typename Type,typename PType>
class CubicIdeal{


public:


// Given a pointer to a cubic order, initialize the ideal as the identity matrix.
CubicIdeal(CubicOrder<Type,PType> * cnfo){
  this->my_order = cnfo;
  for (auto i = 0; i < 3; ++i){
    for (auto j = 0; j < 3; ++j){
      if (i != j){ NTL::clear(this->coeff_matrix[i][j]);}
      else{NTL::set(this->coeff_matrix[i][j]);}
    }
  }
}
/************** Constructor(s) **********************/
CubicIdeal(CubicOrder<Type,PType> * cnfo, const CubicElement<Type,PType> & A,
  const CubicElement<Type,PType> & B, const CubicElement<Type,PType> & C);


/************** Accessors **********************/
CubicOrder<Type, PType> * get_order() const {return my_order;}
const Type get_denom () const {return denom;}
const Type get_coeff(int i, int j) const {
  return coeff_matrix[i][j];
}
void set_order(CubicOrder<Type, PType> * ord) {this->my_order = ord;}

/* @brief Given an integer i from 1..2, assign the input element to be
* the ith element of the Z-basis of the ideal
*/
void copy_basis_element( CubicElement<Type, PType> & elt, int pos) const;
//
inline const CubicElement<Type, PType> * get_gen1()  {
  gen.set_order(my_order);
  gen.assign(coeff_matrix[0][0], coeff_matrix[1][0],coeff_matrix[2][0], this->denom);
  return &gen;
}
inline const CubicElement<Type, PType> * get_gen2()  {
  gen.set_order(my_order);
  gen.assign(coeff_matrix[0][1], coeff_matrix[1][1],coeff_matrix[2][2], this->denom);
  return &gen;}
inline const CubicElement<Type, PType> * get_gen3()  {
  gen.set_order(my_order);
  gen.assign(coeff_matrix[0][2], coeff_matrix[1][2],coeff_matrix[2][2], this->denom);
  return &gen;}


void assign(const CubicElement<Type, PType> g3, const CubicElement<Type, PType> g2, const CubicElement<Type, PType> g1);
void assign(const CubicIdeal<Type, PType> & A);
void assign(const Type U1, const Type X1, const Type Y1, const Type D1,
const Type U2, const Type X2, const Type Y2, const Type D2,
const Type U3, const Type X3, const Type Y3, const Type D3);


/* @brief Changes this ideal into the the identity, which corresponds to O_k itself
*
*/
void make_identity(){
  for (int i =0; i <3; ++i){
    for (int j = 0; j < 3; ++j){
        if (j!= i){
          NTL::clear(coeff_matrix[i][j]);
        }else{
          NTL::set(coeff_matrix[i][j]);
        }
    }
  }
  NTL::set(denom);
}
/*
* @brief returns the numerical value of the ith generator,
* calculated using the 'main' root of the associated CubicOrder
*/
PType get_gen_numeric(int i);


/**
* @brief This returns the puncture of one of generators
* IMPORTANT: This function entirely ignores the denominator of the ideal and
* treats the 3x3 matrix as an ideal I. Further, it implicitly deals with the
* corresponding 1-lattice of I
*/
void puncture(int i,PType & xi, PType & eta);

void puncture_lattice(PType plat[2][2]);
/*
* @brief take coefficient matrix and turns it into an upper triangular matrix.
* Pg. 66 of CFG indicates that the [0][0]th entry is the unique minimal integer
* in the ideal a. Note this does NOT convert to HNF
*/
void make_triangular();

/*
* @brief Returns a bool indicating whether this is an interal ideal or not
* Since the ideal should be in lowest terms, it just sees if denom = 1
*/
bool is_integral(){
  return IsOne(this->denom);
};

bool is_one(){
  if (this->denom != to<Type>(1)){
    return false;
  }

  for (int i =0; i <3; ++i){
    for (int j = 0; j < 3; ++j){
        if ((j!= i) && (coeff_matrix[i][j] != to<Type>(0)) ){
            return false;
        }
        if ((j == i) && (coeff_matrix[i][j] != to<Type>(1)) ){
            return false;
          }
      }
    }
  return true;
}



/**
* @brief function to check equality of ideals (adapted from 1-lattices)
* I am concerned that this check only works for primitive ideals
*/
bool is_equal(const CubicIdeal<Type, PType> & B);

/**
* @brief function to determine whether the ideal is equivalent to B. Not yet implemented
*/
bool is_equivalent(const CubicIdeal<Type, PType> & B);

/**
* @brief Function which calculates the norm of integral ideals.
* At some point should incorporate fractional ideals in a true norm() function
*/
void integral_norm(Type & val );
/**
* @brief function to determine whether the ideal is principal. Not yet implemented
*/
bool is_principal();

/**
* @brief function to determine whether the ideal is prime. Not yet implemented
*/
bool is_prime();




/**
* @brief function to determine whether the ideal is in canonical form. Not yet implemented
*/
bool is_canonical();

/**
* @brief function to obtain an equivalent reduced ideal.
*/
void reduce(CubicIdeal<Type, PType> & rIdeal, CubicElement<Type, PType> & relatingElement);

/**
* @brief function to change the basis into canonical form. Not yet implemented
*
* Important note: To take into account fractional ideals, the assumed state of the
* ideal is as (1/sigma)*J, where sigma is minimal and J is an integral ideal
* Proceed to operate on L1, the 1-lattice corresponding to I,
* L1 = sigma * L2, where L2 corresponds to J
* In hindsight, maybe this isnt as complicated as all that. It's basically
* Find the HNF of J
*/
void make_canonical();

/**
* @brief function to change the idea to a primitive ideal
*
* Important note: The new ideal is NOT the same. It will however be equivalent
* to the original ideal.
* All this function does is check to see if the coeff_matrix entries have a
* non-trivial gcd and divides out the entries.
*/
void make_primitive();

/**
* @brief This function implements Voronoi's algorithm to obtain an adjacent ideal
* Warning this will change the basis of this CubicIdeal to a Voronoi Basis
* @ params[out] B becomes the adjacent ideal
*/
void adjacent_ideal(CubicIdeal<Type, PType> &B, CubicElement<Type, PType> &adj_min, char axis = 'X');




/**
* @brief function to obtain a prepared basis. Relies on the corresponding function
* in the Voronoi.cpp files
*/
void make_prepared(){

  this->puncture_lattice(this->p_lat);
  this->get_order()->get_voronoi()->make_prepared(*this, this->p_lat);
}

/**
* @brief function to convert the basis to a Voronoi basis
*/
void make_voronoi_basis(char axis = 'X');


/*
* @brief This function divides the ideal by the second basis element.
* Generally not useful,It is intended
* for dividing an ideal in voronoi basis form by the adjacent minimum to 1,
* @param[in] B is a CubicIdeal to hold the new lattice, adj_min is a CubicElement
* which will hold the 2nd basis element.
* @param[out] B is the lattice obtained after division
* @pre It is assumed that this CubicIdeal is already in Voronoi form
*/
void divide_adjacent(CubicIdeal<Type, PType> &B, CubicElement<Type, PType> &adj_min);


/*
* @brief Returns a string comprised of the coefficients of the ideal and the denominator
*
*/
std::string toString(){
  string s = "";
  s += (zToString(this->coeff_matrix[0][0]) + " " + zToString(this->coeff_matrix[0][1]) + " " + zToString(this->coeff_matrix[0][2]) + "\n");
  s += (zToString(this->coeff_matrix[1][0]) + " " + zToString(this->coeff_matrix[1][1]) + " " + zToString(this->coeff_matrix[1][2]) + "\n");
  s += (zToString(this->coeff_matrix[2][0]) + " " + zToString(this->coeff_matrix[2][1]) + " " + zToString(this->coeff_matrix[2][2]) + "\n");
  s += "------ \n" + zToString(this->denom);
  return s;
}

/******************** friends *********************/
//friend void mul<Type, PType > (CubicIdeal<Type, PType> & A, const CubicIdeal<Type, PType> & B, const CubicIdeal<Type, PType> & C);
template <typename T, typename PT>
friend class MultiplyStrategyWilliams;
template<typename T, typename PT>
friend class VoronoiMethods;
template<typename T, typename PT>
friend class VoronoiReal;
template<typename T, typename PT>
friend class VoronoiComplex;

template <typename T, typename PT>
friend bool is_equal(CubicIdeal <T,PT> & A, CubicIdeal <T,PT> & B);

static CubicElement<Type, PType> spare_ideal_element;

protected:



// holds a reference to the order in which the ideal sits in
CubicOrder<Type, PType> * my_order;

// a Z basis for the ideal, should be of size 3
// The representation will always be as 3 integral elements all over a common denominator
static Type ci_temp, ci_temp2, ci_temp3;
static Type temp_mat[3][3];
static PType p_temp;
static QQ<Type> rational_temp;

CubicElement<Type, PType> gen;

Type coeff_matrix[3][3];
Type denom{1};
PType p_lat[2][2];

char reduced{2};

/**
* @ brief, returns an indicator whether the ideal is reduced (= 1), not reduced ( = 0),
* or if unknown (= 2)
*/
char is_reduced() const{
  return reduced;
}
void set_reduction_state(char r){
  if (0 <= r && r <= 2){
    reduced = r;
  }
}

void set_coeffs(const CubicElement<Type,PType> & A,
  const CubicElement<Type,PType> & B, const CubicElement<Type,PType> & C);

void normalize();

/* @brief This is a partial validity check to ensure coefficients are assigned
* which yield a full rank Z-module
*
*/
void check_rank();

/** this is used for very specific situations in the prepared basis algorithm.
* No checks are made as to whether this changes any of the ideal's properties
* or violates any representation assumptions.
* Use with discretion.
*/
void flip_basis_element(int i){
  coeff_matrix[0][i] = - coeff_matrix[0][i];
  coeff_matrix[1][i] = - coeff_matrix[1][i];
  coeff_matrix[2][i] = - coeff_matrix[2][i];
}

// **************************************************************************  /
// ******************* Friend classes and functions *************************  /


template <typename T, typename PT>
friend void mul(CubicIdeal<T,PT> & A, const CubicIdeal<T,PT> & B, const CubicIdeal<T,PT> & C);



private:




}; // close class definition


//define statics
template<typename T, typename PT> T CubicIdeal<T,PT>::ci_temp;
template<typename T, typename PT> T CubicIdeal<T,PT>::ci_temp2;
template<typename T, typename PT> T CubicIdeal<T,PT>::ci_temp3;
template<typename T, typename PT> T CubicIdeal<T,PT>::temp_mat[3][3];
template<typename T, typename PT> PT CubicIdeal<T,PT>::p_temp;
template<typename T, typename PT> QQ<T> CubicIdeal<T,PT>::rational_temp;
template<typename T, typename PT> CubicElement<T, PT> CubicIdeal<T,PT>::spare_ideal_element;
#include "../../../src/Cubic/CubicIdeal.cpp"
#endif
