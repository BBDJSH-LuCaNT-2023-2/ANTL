#ifndef ANTL_CUBIC_IDEAL_NF_CPP
#define ANTL_CUBIC_IDEAL_NF_CPP

#include "../../include/ANTL/Cubic/CubicIdealNF.hpp"


template <typename Type, typename PType>
CubicIdealNF<Type,PType> :: CubicIdealNF(const CubicOrderNF<Type,PType> * cnfo, const CubicElement<Type,PType> & A,
  const CubicElement<Type,PType> & B, const CubicElement<Type,PType> & C) : CubicIdeal<Type,PType>::CubicIdeal(cnfo),
  gen1(cnfo, Type(1),Type(0),Type(0),Type(1)),
  gen2(cnfo, Type(1),Type(0),Type(0),Type(1)),
  gen3(cnfo, Type(1),Type(0),Type(0),Type(1)) {
    gen1.assign(A);
    gen2.assign(B);
    gen3.assign(C);
    normalize();
  }

template<typename Type,typename PType>
void CubicIdealNF<Type,PType> :: assign(const CubicElement<Type, PType> g1, const CubicElement<Type, PType> g2, const CubicElement<Type, PType> g3){
  gen1.assign(g1);
  gen2.assign(g2);
  gen3.assign(g3);
  normalize();
}

template<typename Type,typename PType>
void CubicIdealNF<Type,PType> :: assign(const Type U1, const Type X1, const Type Y1, const Type D1,
const Type U2, const Type X2, const Type Y2, const Type D2,
const Type U3, const Type X3, const Type Y3, const Type D3){
  gen1.assign(U1, X1, Y1, D1);
  gen2.assign(U2, X2, Y2, D2);
  gen3.assign(U3, X3, Y3, D3);
  normalize();
}


template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: normalize(){

  //set the new denominator
  mul(this->ci_temp, gen1.get_denom(), gen2.get_denom());
  mul(this->ci_temp, this->ci_temp, gen3.get_denom()); //ci_temp = d1*d2*d3
  mul(this->denom, this->denom, this->ci_temp);

  // Note we should be careful here since mul automatically normalizes after
  // completing the operation. By multiplying by d1*d2*d3, the denominators
  // of each generator should be 1. so there should be no problems normalizing
  // ideal
  mul(gen1, gen1, this->ci_temp);
  mul(gen2, gen2, this->ci_temp);
  mul(gen3, gen3, this->ci_temp);
  this->ci_temp = GCD(GCD(gen1.get_u(), gen2.get_u()),gen3.get_u() );

  if (this->ci_temp != Type(1)){
    this->ci_temp = GCD(GCD(GCD(GCD(GCD(GCD(GCD(
    this->ci_temp,
    gen1.get_x() ),
    gen2.get_x() ),
    gen3.get_x() ),
    gen1.get_y() ),
    gen2.get_y() ),
    gen3.get_y() ),
    this->denom  );

    div(gen1, gen1, this->ci_temp);
    div(gen2, gen2, this->ci_temp);
    div(gen3, gen3, this->ci_temp);
    div(this->denom, this->denom, this->ci_temp);

  }

}


/*
template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_equal(const CubicIdeal<Type, PType> & B){

  std::cout << "... Comparing Lattices ... "<< std::endl;
  Type L1[3][3];
  Type L2[3][3];

  L1[0][0] = this->gen1.get_u();
  L1[1][0] = this->gen1.get_x();
  L1[2][0] = this->gen1.get_u();
  L1[0][1] = this->gen2.get_u();
  L1[1][1] = this->gen2.get_x();
  L1[2][1] = this->gen2.get_y();
  L1[0][2] = this->gen3.get_u();
  L1[1][2] = this->gen3.get_x();
  L1[2][2] = this->gen3.get_y();

  L2[0][0] = B.get_gen1()->get_u();
  L2[1][0] = B.get_gen1()->get_x();
  L2[2][0] = B.get_gen1()->get_y();
  L2[0][1] = B.get_gen2()->get_u();
  L2[1][1] = B.get_gen2()->get_x();
  L2[0][2] = B.get_gen2()->get_y();
  L2[2][1] = B.get_gen3()->get_u();
  L2[1][2] = B.get_gen3()->get_x();
  L2[2][2] = B.get_gen3()->get_y();

  Type denom1 = this.denom;
  Type denom2 = B.get_denom();

  becomeCanonical(L1);
  becomeCanonical(L2);

    this->ci_temp2 = abs(L1[1][0] * L1[2][1] - L1[1][1] * L1[2][0]);
    this->ci_temp3 = abs(L2[1][0] * L2[2][1] - L2[1][1] * L2[2][0]);

    std::cout << "Lattice Invariants: " << denom1<<  " " <<L1[1][0] << " " << (L1[2][1]) << " "
              << L1[1][1] << " " << L1[0][0] << " " << L1[0][1] <<
              "\n" <<
      (L1[1][1]%L1[1][0]) << (L1[0][0]%denom1) <<
      (L1[0][1] - L2[0][1])%denom1<< std::endl;

      std::cout << "Lattice Invariants: " << denom2 <<  " " <<L2[1][0] << " " << (L2[2][1]) << " "
                << L2[1][1] << " " << L2[0][0] << " " << L2[0][1] << "\n" <<
      (L2[1][1]%L1[1][0]) << (L2[0][0]%denom1) <<
      ( ((L1[1][1] - L2[1][1]) / L1[1][0] )*L1[0][0]) % denom1<< std::endl;


    if (   ( denom1 != denom2   )     // matching denominators
        || ( this->ci_temp2 != this->ci_temp3 )               // e values match
        || (L1[2][1] != L2[2][1]) ){  // g values match

        std::cout << "first equality check fails" << std::endl;
        std::cout << "----------------"<< std::endl;
        return false;
    }//close if clause
    else if (
         (  ( (L1[1][1] - L2[1][1]) % L1[1][0] ) != 0  )

      || (  ( (L1[0][0] - L2[0][0]) %denom1) != 0 )
      || (   ( (L1[0][1] - L2[0][1])%denom1)
          != ( ( (L1[1][1] - L2[1][1])*L1[0][0] / L1[1][0] ) % denom1)   )
      ) //close else if clause
    {

        std::cout << "second equality check fails" << std::endl;

        std::cout << "----------------"<< std::endl;
        return false;
    }
    else{
      return true;
    }
} //close function is_equal
*/




template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_equivalent(const CubicIdeal<Type, PType> & B){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
Type CubicIdealNF<Type, PType> :: norm(){
  std::cout <<"undefined" << std::endl;
  return Type(-1);
}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_principal(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_prime(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
bool CubicIdealNF<Type, PType> :: is_canonical(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: reduce(){
  std::cout <<"undefined" << std::endl;
}








#endif
