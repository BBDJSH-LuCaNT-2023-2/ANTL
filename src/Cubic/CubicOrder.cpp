#ifndef ANTL_CUBIC_ORDER_CPP
#define ANTL_CUBIC_ORDER_CPP

#include "../../include/ANTL/Cubic/CubicOrder.hpp"
#include <algorithm>

// constructor definitions
template <typename Type, typename PType>
CubicOrder<Type, PType> :: CubicOrder(polynomial<Type> const &poly, std::string s) {

    this->discriminant = discriminant_bcf(poly);

    if (s.compare("Williams") == 0){
        if (!(this->williams)){
          this->williams.reset();
          this->williams = std::make_shared<MultiplyStrategyWilliams<Type, PType>>();
        }

        this->m_method = std::static_pointer_cast< IdealMultiplicationStrategy<Type, PType>>(this->williams);
    }
    else{
      std::cout << "Invalid multiplication strategy provided" << std::endl;
    }
    defining_IBCF = poly;

    // this initiate the fund unit strategy, we only have one for now, but the if statement will
    // allow us to switch strategies possibly or just use the setter method
    if (true){
        if (!(this->voronoi_basic)){
          this->voronoi_basic.reset();
          this->voronoi_basic = std::make_shared<BasicVoronoi<Type, PType>>();
        }

        this->unit_strat = std::static_pointer_cast< FundUnitStrategy<Type, PType>>(this->voronoi_basic);
    }


    if (this->discriminant < 0){
      //this->vmethods = std::static_pointer_cast< VoronoiMethods<Type, PType>>( std::make_shared<VoronoiComplex<Type,PType> >() );
      this->vcomplex.reset();
      this->vcomplex = std::make_shared<VoronoiComplex<Type, PType> >();
      this->vmethods = std::static_pointer_cast <VoronoiMethods<Type, PType> > (this->vcomplex);
    }
    else{
      this->vreal.reset();
      this->vreal = std::make_shared<VoronoiReal<Type, PType> >();
      this->vmethods = std::static_pointer_cast <VoronoiMethods<Type, PType> > (this->vreal);
    }

    set_roots();
    set_integral_basis();
    set_mul_table();


}


template <typename Type, typename PType>
bool CubicOrder<Type, PType> :: is_equal(const CubicOrder<Type, PType> &CO2) const{

  return (this->defining_IBCF[0] == CO2.get_IBCF()[0] && this->defining_IBCF[1] == CO2.get_IBCF()[1] && this->defining_IBCF[2] == CO2.get_IBCF()[2] && this->defining_IBCF[3] == CO2.get_IBCF()[3]);
}



// definitions for the getters which are not automatically computed
template <typename Type, typename PType>
ZZ CubicOrder<Type, PType> :: get_class_number(){
  if (!IsZero(class_number)){
    return class_number;
  } else{
    set_class_number();
    return class_number();
  }
}


template <typename Type, typename PType>
PType CubicOrder<Type, PType> :: get_regulator(){
  if (regulator != 0){
    return regulator;
  } else{
    set_regulator();
    return regulator();
  }
}


template <typename Type, typename PType>
std::vector<Type> CubicOrder<Type, PType> :: get_class_group(){
  if (cg_structure.size() != 0){
    return cg_structure;
  } else{
    set_class_group();
    return cg_structure();
  }
}

template <typename Type, typename PType>
void CubicOrder<Type, PType> ::  roots_swap_position(int p1, int p2){
  if (this->discriminant > 0) {
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
void CubicOrder<Type, PType> :: get_real_value(PType & newVal, Type & U, Type & X, Type & Y, Type & D, int conj) {
    std::cout << "using CubicOrders real function" << std::endl;
    NTL::mul(newVal, to<PType>(X), this->get_rho1());
    NTL::add(newVal, newVal, to<PType>(U));
    NTL::mul(this->order_temp, to<PType>(Y), this->get_rho2());
    NTL::add(newVal, newVal, this->order_temp);
    NTL::div(newVal, newVal, to<PType>(D));


}



// *********************Protected member methods  *************************** //
template <typename Type, typename PType>
void CubicOrder<Type, PType> :: set_roots(){
  int type = cardano<Type, PType>(this->defining_IBCF ,this->root_list );
  PType x1;
  abs(x1, this->root_list[0]);
  if ( (this->discriminant > 0) && (x1 - to<PType>(1.0) <= to<PType>(0)) ){
    //roots_swap_position(this->root_list[0], this->root_list[1]);
    roots_swap_position(0, 1);
  }
}


template <typename Type, typename PType>
void CubicOrder<Type, PType> :: set_mul_table( ){

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
void CubicOrder<Type, PType> :: set_integral_basis(){

  // rho1 = a*delta
  NTL::mul(this->rho1, this->root_list[0], to<PType>(this->defining_IBCF[3]) );

  // rho2 = a*delta + b
  add(this->rho2, this->rho1, to<PType>(this->defining_IBCF[2]) );
  // rho2 = a*delta*delta + b*delta
  NTL::mul(this->rho2, this->rho2, this->root_list[0]);

}


template <typename Type, typename PType>
void CubicOrder<Type, PType> :: set_class_number( ){

    std::cout << "set class number"<< std::endl;
}


template <typename Type, typename PType>
void CubicOrder<Type, PType> :: set_class_group( ){

    std::cout << "set class group"<< std::endl;
}

/*
template <typename Type, typename PType>
void CubicOrder<Type, PType> :: compute_fundamental_unit(){

  // initialize matrix L0 as the identity

  CubicIdeal<Type, PType> L_identity = CubicIdeal<Type, PType>(this);
  CubicIdeal<Type, PType> L1 = CubicIdeal<Type, PType>(this);
  CubicIdeal<Type, PType> L2 = CubicIdeal<Type, PType>(this);

  CubicElement<Type, PType> epsilon = CubicElement<Type, PType>(this);
  CubicElement<Type, PType> adj_minimum = CubicElement<Type, PType>(this);
  int rounds=0;
  bool complete_cycle=false;

  do {

    //////////////////////////////////////
    // For testing purposes only: creates test-lattice Ltest
    // copy L to Ltest
    //for (int i = 0; i < 3; ++i){
    //  for (int j = 0; j < 2; ++j){
    //    Ltest.coefficientMatrix[i][j] = L.coefficientMatrix[i][j];
    //  }
    //}
    //Ltest.mainDenominator = L.mainDenominator;
    /////////////////////////////////////

    std::cout << " *************************************** " << std::endl;
    std::cout << " *************************************** " << std::endl;

        std::cout << "Fund.Unit Iteration Number --- " << rounds << std::endl;
        ++rounds; //prints and increases counter.


        // compute the VoronoiBasis (1, theta_g, theta_h),
        //  where theta_g is the minima adjacent to 1 in L
        //L.make_voronoi_basis();

//std::cout << "FundamentalUnit2: After voronoiBasis: " << std::endl;
//std::cout << L.coefficientMatrix[0][0] << "   " << L.coefficientMatrix[0][1] << " " << std::endl;
//std::cout << L.coefficientMatrix[1][0] << "   " << L.coefficientMatrix[1][1] << " " << std::endl;
//std::cout << L.coefficientMatrix[2][0] << "   " << L.coefficientMatrix[2][1] << " " << std::endl;
//std::cout << L.mainDenominator << std::endl;
//std::cout << " ****************************** " << std::endl;
//std::cout << "FundamentalUnit2: Is VB the same? "<< std::endl;
//        if (compareLattice3(Ltest,L)!= 1){
//          std::cout << "ERROR: VoronoiBasis has returned a different lattice. ABORT";
//          break;
//        }
//        else{
//          std::cout << "VoronoiBasis has returned the same lattice, continuing:";
//        }
//std::cout << " ****************************** " << std::endl;


        // If proper, after this call, L2 is the adjacent lattice, and L1. spare_ideal_element
        // will contain the relative minima adjacent to 1.

        L1.adjacent_ideal(L2, adj_minimum );

        // mulitiply epsilon by the min adjacent to 1.
        ::mul(epsilon, epsilon, adj_minimum );

        std::swap(L1, L2);

        complete_cycle = ::is_equal(L1, L_identity);

    } while( !complete_cycle );

    std::cout << "Fundamental Unit computed: "<< std::endl;
  this->fundamentalUnits.push_back(epsilon);
}
*/

template <typename Type, typename PType>
void CubicOrder<Type, PType> :: set_regulator( ){
    std::cout << "set_regulator not yet implemented"<< std::endl;

}


// ***********************Friend function definitions
template <typename T, typename PT>
bool is_equal(const CubicOrder<T, PT> &CO1, const CubicOrder<T, PT> &CO2) {
  return (
    CO1.defining_IBCF[0] == CO2.defining_IBCF[0] && CO1.defining_IBCF[1] == CO2.defining_IBCF[1] && CO1.defining_IBCF[2] == CO2.defining_IBCF[2] && CO1.defining_IBCF[3] == CO2.defining_IBCF[3]
  );
}


#endif
