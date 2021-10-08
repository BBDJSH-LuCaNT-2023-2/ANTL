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
    set_unit_rank();
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
    return regulator;
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




//template <typename Type, typename PType>
//void CubicOrder<Type, PType> :: get_real_value(PType & newVal, Type & U, Type & X, Type & Y, Type & D, int conj) {
//    std::cout << "using CubicOrders real function" << std::endl;
//    NTL::mul(newVal, to<PType>(X), this->get_rho1());
//    NTL::add(newVal, newVal, to<PType>(U));
//    NTL::mul(this->order_temp, to<PType>(Y), this->get_rho2());
//    NTL::add(newVal, newVal, this->order_temp);
//    NTL::div(newVal, newVal, to<PType>(D));
//}



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

template <typename Type, typename PType>
void CubicOrder<Type, PType> :: set_unit_rank( ){

    if(this->discriminant < 0){
      this->unit_rank = 1;
    }else{
      this->unit_rank = 2;
    }
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
void CubicOrder<Type, PType> :: standard_form(Type & E, Type& G){


  Type temp_var1, temp_var2;
  mul(E, this->get_coeff(2),this->get_coeff(2));    // b^2
  mul(temp_var1, Type(3), this->get_coeff(3));
  mul(temp_var1, temp_var1, this->get_coeff(1));    // 3ac
  sub (E, E, temp_var1);                            // b^2 - 3ac   // A

  mul(G, E, this->get_coeff(2));                    // Ab

  mul(E, E, Type(-3));                              // E = -3A

  mul(G, G, Type(2));                               // G = 2Ab

  mul(temp_var1, this->get_coeff(1), this->get_coeff(2)); // bc
  mul(temp_var2, this->get_coeff(0), this->get_coeff(3)); // ad
  mul(temp_var2, Type(9), temp_var2);                      // 9ad
  sub(temp_var1, temp_var1, temp_var2);                   // temp_var1 = B
  mul(temp_var1, this->get_coeff(3), temp_var1);          // temp_var1 = aB
  mul(temp_var1, temp_var1, Type(3));                     // 3aB
  sub(G, G, temp_var1);                                   // G = 2Ab - 3aB

}

// We have to denote the
template <typename Type, typename PType>
int CubicOrder<Type, PType> :: splitting_type(Type p){

  // make sure that the standard form has been computed before
  if ( (this->E == 0) && (this->G == 0) ){
    this->standard_form(E,G);
  }
  Type E_t = this->E;
  Type G_t = this->G;
  int sp = 0; // should hold the valuation of p in discriminant of standard form f, v_p(D_f)
  Type Dp; // holds disc(f) / p^sp

  // Compute the discriminant of standard form poly
  sqr(Dp, E_t);
  mul(Dp, Dp, E_t);
  mul(Dp, Dp, Type(4));         // 4E^3
  sqr(sp, G_t);
  mul(sp, sp, Type(27));        // 27G^2
  sub(Dp, Dp, sp);              // 4E^3-27G^2
  while(Dp % p == Type(0)){
    sp+= 1;
    div(Dp, Dp, p);
  }
  //  0 = p^3
  //  1 = p * q^2
  //  2 = p*q
  //  3 = p* q * r
  //  4 = p
  if (p == Type(2)){
    ///////////////////////////////
    if( (E_t%2 ==0 ) && (G_t%2 == 0) ){
      do{
        div(E_t, E_t, p);
        div(G_t,G_t, p);
        if( (E_t% p != Type(0)) && (G_t %p == Type(0)) ){
          return 1;
        }
        if( ((E_t% p == Type(0)) && (G_t %p != Type(0))) || ( (E_t% p != Type(0)) && (G_t %p != Type(0)) ) ){
          return 0;
        }
      }while(true);
    }
    ///////////////////////////////
    else if ((E_t%2 ==0 ) && (G_t%2 != 0)){
      return 2;
    }
    ///////////////////////////////
    else if( (E_t%2 !=0 ) && (G_t%2 == 0) ){
      if (sp%2 == 1){return 1;}
      else {
        if (Dp % 4 == 3){return 1;}
        if (Dp % 8 == 5){return 2;}
        if (Dp % 8 == 1){return 3;}
      }
    }
    ///////////////////////////////
    else {
      // this is the case where 2 divides neither G nor E_t
      return 4;
    }
  } // end of case p = 2
  //////////////////////////////////////////////////////////////
  //                          p = 3                           //
  //////////////////////////////////////////////////////////////
  else if ( p == Type(3)){
    ///////////////////////////////
    if( (E_t%p ==0 ) && (G_t%p == 0) ){
      do{
        div(E_t, E_t, p);
        div(G_t,G_t, p);
        if( (E_t% p != Type(0)) && (G_t %p == Type(0)) ){
          return 1;
        }
        else if( ((E_t% p == Type(0)) && (G_t %p != Type(0))) || ( (E_t% p != Type(0)) && (G_t %p != Type(0)) ) ){
          return 0;
        }
      }while(true);
    }
    else if((E_t%p ==0 ) && (G_t%p != 0)){
        sqr(G_t, G_t);
      if (E_t%9 != 3){
        add(E_t, E_t, Type(1));
        sub(G_t, G_t,E_t);
        if (G_t % 9 == 0){
          return 1;
        }else{
          return 4;
        }
      }
      else {
        add(E_t, E_t, Type(1));
        sub(G_t, G_t,E_t);
        if (G_t % 27 == 0){
          if (sp%2 == 1){return 1;}
          else if (Dp %3 == 2){return 2;}
          else if ( (sp == 6) && (Dp%3 ==1 ) ){return 4;}
          else{return 3;}
        }
        else {
          return 0;
        }
      }
    }
    ///////////////////////////////
    else {
      if (E_t%p == Type(2)){
        return 2;
      }
      else if (G_t % p != 0){
        return 4;
      }
      else{
        return 3;
      }

    }
  } // end p = 3 case
  //////////////////////////////////////////////////////////////
  //                          p >= 5                          //
  //////////////////////////////////////////////////////////////
  else{
    // p divides both E and G
    if( (E_t%p ==0 ) && (G_t%p == 0) ){
      do{
        div(E_t, E_t, p);
        div(G_t,G_t, p);

        // val(G) > val(E)
        if( (E_t% p != Type(0)) && (G_t %p == Type(0)) ){
          return 1;
        }
        // Case when val(G) <= val(E)
        else if( ((E_t% p == Type(0)) && (G_t %p != Type(0))) || ( (E_t% p != Type(0)) && (G_t %p != Type(0)) ) ){
          return 0;                                           //
        }
      }while(true);
    }
    // p divides E but not G
    else if((E_t%p ==0 ) && (G_t%p != 0)){
        if (p%3 == 2){return 2;}
        else{                   // https://math.stackexchange.com/questions/127251/when-is-a-not-a-cube-mod-p
          Type legendre, legexp;
          sub(legexp, p, Type(1));
          div(legexp, p, Type(3));        // p-1/3
          NTL::PowerMod(legendre, G_t,legexp, p);
          if ( (IsOne(legendre)) )
            {return 3;}
          else
            {return 4;}
        }
    }
    //Exponentiation<Type> bexp;
    // p divides G but not E

    else if((E_t%p !=0 ) && (G_t%p == 0)){
      Type legendre, legexp;
      sub(legexp, p, Type(1));
      div(legexp, legexp, 2);
        NTL::PowerMod(legendre, G_t,legexp, p);
        if ((legendre == 1)){return 3;}
        else{return 2;}
    }
    // p does not divide G or E
    else{
      Type legendre, legexp;
      sub(legexp, p, Type(1));
      div(legexp, legexp, 2);                       // legexp = (p-1)/2
      if (sp %2 == 1){return 1;}                    // check sp odd
      else{
        NTL::PowerMod(legendre, Dp,legexp, p);
        if (legendre == 1){
          bool reduc;
          Type result;
          int i = 0;
          while (i < p && !(IsZero(result)) ){
            eval_cubic_mod_p(result, Type(i), this->defining_IBCF[0],this->defining_IBCF[1],this->defining_IBCF[2], this->defining_IBCF[3], p);
            i++;
          }
          if (IsZero(result)){return 3;}   /// just plug in all values? Better way?
          else {return 4; }
        }
        else {return 2; }
      }
    }
  } // end p > 3 case
} // end splitting_type function


// ***********************Friend function definitions
template <typename T, typename PT>
bool is_equal(const CubicOrder<T, PT> &CO1, const CubicOrder<T, PT> &CO2) {
  return (
    CO1.defining_IBCF[0] == CO2.defining_IBCF[0] && CO1.defining_IBCF[1] == CO2.defining_IBCF[1] && CO1.defining_IBCF[2] == CO2.defining_IBCF[2] && CO1.defining_IBCF[3] == CO2.defining_IBCF[3]
  );
}


#endif
