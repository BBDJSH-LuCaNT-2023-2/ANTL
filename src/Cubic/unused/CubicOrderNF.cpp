#ifndef ANTL_CUBIC_ORDERNF_CPP
#define ANTL_CUBIC_ORDERNF_CPP

#include "../../include/ANTL/Cubic/CubicOrderNF.hpp"



// *********************** Public method definitions

template <typename Type, typename PType>
CubicOrderNF<Type, PType> :: CubicOrderNF( polynomial<Type> const &poly)
  : CubicOrder<Type,PType>::CubicOrder(poly)  {

    this->discriminant = discriminant_bcf(poly);

    set_roots();
    set_integral_basis();
    set_mul_table();

}


template <typename Type, typename PType>
bool CubicOrderNF<Type, PType> :: is_equal(const CubicOrder<Type, PType> &CO2) const{

  return (this->defining_IBCF[0] == CO2.get_IBCF()[0] && this->defining_IBCF[1] == CO2.get_IBCF()[1] && this->defining_IBCF[2] == CO2.get_IBCF()[2] && this->defining_IBCF[3] == CO2.get_IBCF()[3]);
}


template <typename Type, typename PType>
void CubicOrderNF<Type, PType> ::  roots_swap_position(int p1, int p2){
  if (this->discriminant < 0) {
    // swap roots and recalculate the integral basis.
    std::swap(this->root_list[p1], this->root_list[p2]);
    set_integral_basis();
  }
  else{
    std::cout << "Attempting to swap root order for a complex cubic order: root_list is a fixed order: Real root, real part, imaginary part. \
    no action taken." << std::endl;
  }
};

template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: mul(CubicIdeal<Type,PType> & A, const CubicIdeal<Type,PType> & B, const CubicIdeal<Type,PType> & C){
  if (!(A.get_order()->is_equal( (*B.get_order()) ) ) || (!(A.get_order()->is_equal( (*this) ) ) )){
    std::cout << "Order mismatch in ideal mul" << std::endl;
  }
  (this->mul_method)->multiply(A,B,C);
}




// *********************Protected member methods
template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_roots(){
  int type = SolveP3<PType>(this->root_list, PType(this->defining_IBCF[2]/this->defining_IBCF[3]),
    PType(this->defining_IBCF[1]/this->defining_IBCF[3]), PType(this->defining_IBCF[0]/this->defining_IBCF[3]) );

  if ( (this->discriminant > 0) &&(ANTL::abs(this->root_list[0])-PType(1.0) <= PType(0))){
    roots_swap_position(this->root_list[0], this->root_list[1]);
  }
}


template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_mul_table( ){

  // fill col1, corresponds to rho1^2
  mul_table[0][0] = 0;                                          // 0
  mul_table[1][0] = -(this->defining_IBCF[2]);                  // -b
  mul_table[2][0] = this->defining_IBCF[3];                     // a

  //col2 corresponds to rho1*rho2
  NTL::mul(mul_table[0][1], -this->defining_IBCF[3],this->defining_IBCF[0]);    //-ad
  //mul_table[0][1] = - (this->defining_IBCF[0])*(this->defining_IBCF[3]);
  mul_table[1][1] = -this->defining_IBCF[1];                                    // -c
  mul_table[2][1] = 0;
  NTL::mul(mul_table[0][2], -this->defining_IBCF[2],this->defining_IBCF[0] );   // -bd
  mul_table[1][2] = -this->defining_IBCF[0] ;                                   // -d
  mul_table[2][2] = -this->defining_IBCF[1] ;                                   // -c
  //std::cout << mul_table[2][0] << mul_table[2][1] << mul_table[2][2]<< std::endl;
}


template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_integral_basis(){

  // rho1 = a*delta
  NTL::mul(this->rho1, this->root_list[0], PType(this->defining_IBCF[3]) );

  // rho2 = a*delta + b
  add(this->rho2, this->rho1, this->defining_IBCF[2]);
  // rho2 = a*delta*delta + b*delta
  NTL::mul(this->rho2, this->rho2, this->root_list[0]);

}


template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_class_number( ){

    std::cout << "set class number"<< std::endl;
}


template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_class_group( ){

    std::cout << "set class group"<< std::endl;
}


template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_regulator( ){
    std::cout << "hi"<< std::endl;

}


// ***********************Friend function definitions
template <typename T, typename PT>
bool is_equal(const CubicOrderNF<T, PT> &CO1, const CubicOrderNF<T, PT> &CO2) {
  return (
    CO1.defining_IBCF[0] == CO2.defining_IBCF[0] && CO1.defining_IBCF[1] == CO2.defining_IBCF[1] && CO1.defining_IBCF[2] == CO2.defining_IBCF[2] && CO1.defining_IBCF[3] == CO2.defining_IBCF[3]
  );
}




#endif
