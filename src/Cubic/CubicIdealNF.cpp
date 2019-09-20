#ifndef ANTL_CUBIC_IDEAL_NF_CPP
#define ANTL_CUBIC_IDEAL_NF_CPP

#include "../../include/ANTL/Cubic/CubicIdealNF.hpp"


template <typename Type, typename PType>
CubicIdealNF<Type,PType> :: CubicIdealNF(const CubicOrderNF<Type,PType> * cnfo, const CubicElementNF<Type,PType> & A,
  const CubicElementNF<Type,PType> & B, const CubicElementNF<Type,PType> & C) : CubicIdeal<Type,PType>::CubicIdeal(cnfo),
  gen1(cnfo, Type(1),Type(0),Type(0),Type(1)),
  gen2(cnfo, Type(1),Type(0),Type(0),Type(1)),
  gen3(cnfo, Type(1),Type(0),Type(0),Type(1)) {
    gen1.assign(A);
    gen2.assign(B);
    gen3.assign(C);
    normalize();
  }

template<typename Type,typename PType>
void CubicIdealNF<Type,PType> :: assign(const CubicElementNF<Type, PType> g1, const CubicElementNF<Type, PType> g2, const CubicElementNF<Type, PType> g3){
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

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: become_canonical(){
  std::cout <<"undefined" << std::endl;


    Type g, r, s;          // store gcd g and integer multipliers r,s which
                                 // satisfy r*m_32 + s*m_33 = g

    Type dummy, j,k,       //dummy is extraneous, j,k satisfy jr-ks = 1
                q,                // transitory value
                minorDet;         // determinant of minor matrix

    Type canonMatrix[3][3]; //variable for the canonical matrix

    // updateLattice removes common factors in denominator and matrix entries
    // Commented out since CubicIdeals should always be in lowest terms
    // updateLattice(B);

    // Checks whether the properties of canonical basis are already satisfied
      if(
        ( this->gen2.get_y() == 0 ) &&
        ( 0 <= this->gen2.get_u() && this->gen2.get_u() < this->gen1.get_u()) &&
        ( 0 <= this->gen3.get_u() && this->gen3.get_u() < this->gen1.get_u()) &&
        ( 0 <= this->gen3.get_x() && this->gen3.get_x() < this->gen2.get_x())
      ) //close if clause
      {
        // do nothing if already canonical
      }
      else{
      //computes the determinant of the lower-right 2x2 minor
      minorDet = (this->gen2.get_x() * this->gen3.get_y()
                  - this->gen3.get_x() * this->gen2.get_y());

      //From the NTL Library, EEA
      XGCD(g, r, s ,this->gen2.get_y(), this->gen3.get_y());

      XGCD(dummy, j,k,r,-s);    // jr - ks = 1, dummy is always 1

      q = -(k*this->gen2.get_y() + j*this->gen3.get_y());
      q = q/g;         // computes transitory value

      dummy = k+r*q;  //reusing dummy to save a couple ops


      canonMatrix[0][0] =  (dummy) * this->gen2.get_u() + (j+q*s)*this->gen3.get_u();
      canonMatrix[1][0] =  -minorDet/g;  // denoted  -epsilon/g in cubic book
      canonMatrix[2][0] =  0;
      canonMatrix[0][1] = r * this->gen2.get_u() +  s*this->gen3.get_u();
      canonMatrix[1][1] = r * this->gen2.get_x() +  s*this->gen3.get_x();
      canonMatrix[2][1] = g;

      // if -minorDet is negative, we need to flip the middle row to get e/g, where
      //  e is the absolute value of minorDet
      if (canonMatrix[1][0] < 0){
          canonMatrix[1][0] = -canonMatrix[1][0];
          canonMatrix[0][0] = - canonMatrix[0][0];
      }


      //Ensures that canonMatrix[1][1] (this is m_23 in the cubic book) satisfies
      // 0 <= m_23 < e/g )
      while ( (canonMatrix[1][1] >= canonMatrix[1][0]) || (canonMatrix[1][1] < 0) ){
          if (canonMatrix[1][1] < 0){
              canonMatrix[0][1] += canonMatrix[0][0];
              canonMatrix[1][1] += canonMatrix[1][0];
          }
          else if (canonMatrix[1][1] >= canonMatrix[1][0]){
              canonMatrix[0][1] -= canonMatrix[0][0];
              canonMatrix[1][1] -= canonMatrix[1][0];
          }
      }

      //ensures that m_12 is in the appropriate range
      while ( (canonMatrix[0][0] >= this->gen1.get_u()) || (canonMatrix[0][0] < 0) ){
          if (canonMatrix[0][0] < 0){
              canonMatrix[0][0] += this->gen1.get_u();
          }
          else if (canonMatrix[0][0] >= this->gen1.get_u()){
              canonMatrix[0][0] -= this->gen1.get_u();
          }
      }
      //ensures m_13 is in the appropriate range
      while ( (canonMatrix[0][1] >= this->gen1.get_u()) || (canonMatrix[0][1] < 0) ){
          if (canonMatrix[0][1] < 0){
              canonMatrix[0][1] += this->gen1.get_u();
          }
          else if (canonMatrix[0][1] >= this->gen1.get_u()){
              canonMatrix[0][1] -= this->gen1.get_u();
          }
      }

      gen1.assign(canonMatrix[0][0],canonMatrix[1][0],canonMatrix[2][0], Type(1));
      gen1.assign(canonMatrix[0][1],canonMatrix[1][1],canonMatrix[2][1], Type(1));

      this->normalize();
    }


}






template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: become_prepared(){
  std::cout <<"undefined" << std::endl;
}

template<typename Type,typename PType>
void CubicIdealNF<Type, PType> :: voronoi(){
  std::cout <<"undefined" << std::endl;
}
#endif
