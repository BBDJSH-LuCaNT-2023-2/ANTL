#ifndef ANTL_CUBIC_ORDERNF_CPP
#define ANTL_CUBIC_ORDERNF_CPP

#include "../../include/ANTL/Cubic/CubicOrderNF.hpp"



template <typename Type, typename PType>
CubicOrderNF<Type, PType> :: CubicOrderNF( polynomial<Type> const &poly)
  : CubicOrder<Type,PType>::CubicOrder(poly)  {

    this->discriminant = discriminant_bcf(poly);

    int type = SolveP3<PType>(this->root_list, PType(poly[2]/poly[3]),
      PType(poly[1]/poly[3]), PType(poly[0]/poly[3]) );

    set_integral_basis();
    set_mul_table();
}


template <typename Type, typename PType>
bool CubicOrderNF<Type, PType> :: is_equal(const CubicOrder<Type, PType> &CO2) const{

  return (this->defining_IBCF[0] == CO2.get_IBCF()[0] && this->defining_IBCF[1] == CO2.get_IBCF()[1] && this->defining_IBCF[2] == CO2.get_IBCF()[2] && this->defining_IBCF[3] == CO2.get_IBCF()[3]);
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
  std::cout << mul_table[2][0] << mul_table[2][1] << mul_table[2][2]<< std::endl;
}

template <typename Type, typename PType>
void CubicOrderNF<Type, PType> :: set_integral_basis( ){

  std::cout << "set integral basis"<< std::endl;
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

template <typename T, typename PT>
bool is_equal(const CubicOrder<T, PT> &CO1, const CubicOrder<T, PT> &CO2) {
  return (
    CO1.defining_IBCF[0] == CO2.defining_IBCF[0] && CO1.defining_IBCF[1] == CO2.defining_IBCF[1] && CO1.defining_IBCF[2] == CO2.defining_IBCF[2] && CO1.defining_IBCF[3] == CO2.defining_IBCF[3]
  );
}




#endif
