#ifndef ANTL_CUBIC_IDEAL_CPP
#define ANTL_CUBIC_IDEAL_CPP

#include "../../include/ANTL/Cubic/CubicIdeal.hpp"






template <typename Type, typename PType>
CubicIdeal<Type,PType> :: CubicIdeal(CubicOrder<Type,PType> * cnfo, const CubicElement<Type,PType> & A,
  const CubicElement<Type,PType> & B, const CubicElement<Type,PType> & C):
    gen(cnfo, Type(1),Type(0),Type(0),Type(1)) {
    this->my_order = cnfo;
    set_coeffs(A,B,C);
  }


template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: set_coeffs(const CubicElement<Type,PType> & A,
  const CubicElement<Type,PType> & B, const CubicElement<Type,PType> & C){
    mul(this->ci_temp, B.get_denom(), C.get_denom());
    mul(coeff_matrix[0][0],this->ci_temp, A.get_u()  );
    mul(coeff_matrix[1][0],this->ci_temp, A.get_x() );
    mul(coeff_matrix[2][0],this->ci_temp, A.get_y() );
    mul(this->ci_temp, A.get_denom(), C.get_denom());
    mul(coeff_matrix[0][1],this->ci_temp, B.get_u() );
    mul(coeff_matrix[1][1],this->ci_temp, B.get_x() );
    mul(coeff_matrix[2][1],this->ci_temp, B.get_y() );
    mul(this->ci_temp, B.get_denom(), A.get_denom());
    mul(coeff_matrix[0][2],this->ci_temp, C.get_u() );
    mul(coeff_matrix[1][2],this->ci_temp, C.get_x() );
    mul(coeff_matrix[2][2],this->ci_temp, C.get_y() );

    mul(this->denom, this->ci_temp, C.get_denom());     // d1*d2*d3

    //std::cout << "set coeffs" << std::endl;
    //std::cout << coeff_matrix[0][0] << " " << coeff_matrix[0][1] << " " << coeff_matrix[0][2] << std::endl;
    //std::cout << coeff_matrix[1][0] << " " << coeff_matrix[1][1] << " " << coeff_matrix[1][2] << std::endl;
    //std::cout << coeff_matrix[2][0] << " " << coeff_matrix[2][1] << " " << coeff_matrix[2][2] << std::endl;

    make_canonical();
}


template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: check_rank(){

}



template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: assign(const CubicElement<Type, PType> g1, const CubicElement<Type, PType> g2, const CubicElement<Type, PType> g3){
  set_coeffs(g2,g2,g3);
}

template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: assign(const Type U1, const Type X1, const Type Y1, const Type D1,
const Type U2, const Type X2, const Type Y2, const Type D2,
const Type U3, const Type X3, const Type Y3, const Type D3){


  mul(this->ci_temp, D2, D3);
  mul(coeff_matrix[0][0],this->ci_temp, U1  );
  mul(coeff_matrix[1][0],this->ci_temp, X1 );
  mul(coeff_matrix[2][0],this->ci_temp, Y1 );
  mul(this->ci_temp, D1, D3);
  mul(coeff_matrix[0][1],this->ci_temp, U2 );
  mul(coeff_matrix[1][1],this->ci_temp, X2 );
  mul(coeff_matrix[2][1],this->ci_temp, Y2 );
  mul(this->ci_temp, D2, D1);
  mul(coeff_matrix[0][2],this->ci_temp, U3 );
  mul(coeff_matrix[1][2],this->ci_temp, X3 );
  mul(coeff_matrix[2][2],this->ci_temp, Y3 );

  mul(this->denom, this->ci_temp, D3);     // d1*d2*d3
  this->make_canonical();
}
template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: assign(const CubicIdeal<Type, PType> & A){
  this->set_order(A.get_order());
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j< 3; ++j){
      this->coeff_matrix[i][j] = A.get_coeff(i,j);
    }
  }
  this->denom = A.get_denom();
  this->reduced = A.is_reduced();
}

template<typename Type,typename PType>
PType CubicIdeal<Type,PType> :: get_gen_numeric(int i){
  return this->coeff_matrix[0][i] + this->coeff_matrix[1][i]*(this.get_order()->get_rho1()) + this->coeff_matrix[2][i]*(this.get_order()->get_rho2());
}


// Note I'm ignoring the ideal's true denominator and using
template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: puncture(int i, PType & xi, PType & eta){

  // code to compute xi
  // note that this assigns the denominator of gen to be coeff_matrix[0][0]
  gen.set_order(this->get_order());
  gen.assign(this->coeff_matrix[0][i], this->coeff_matrix[1][i], this->coeff_matrix[2][i], this->coeff_matrix[0][0]);

  gen.get_real_value(xi);

  mul(xi, xi, to<PType>(3));
  gen.trace( this->rational_temp);

  //std::cout << gen.get_u() << " " << gen.get_x() << " " << gen.get_y() << " " << gen.get_denom() << std::endl;
  //std::cout << xi << std::endl;
  //std::cout << this->rational_temp.getNumerator() << " / " << this->rational_temp.getDenominator() << std::endl;
  div(this->p_temp, to<PType>(this->rational_temp.getNumerator()), to<PType>(this->rational_temp.getDenominator()) );
  sub(xi, xi, this->p_temp);
  mul(xi, xi, 0.5);
  // now compute eta
  mul(this->p_temp, to<PType>(this->coeff_matrix[2][i]) ,this->get_order()->get_root1()); // y*delta
  eta = to<PType>(this->coeff_matrix[1][i]);                                      // x
  sub(eta, eta, this->p_temp);                                                // x - y*delta

  if ( this->get_order()->is_complex() ){
    mul(eta, eta, to<PType>(this->get_order()->get_coeff(3)));      // should mult by a    // a*(x - y*delta)
    mul(eta, eta, this->get_order()->get_root3());                              // above * (imaginary part of roots)

  }
  else{
    sub(this->p_temp, this->get_order()->get_root2(), this->get_order()->get_root3());
    mul(eta, eta, this->p_temp );
    mul(eta, eta, 0.5);

  }
  div(eta, eta, to<PType>(this->coeff_matrix[0][0]) );

};

template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: puncture_lattice(PType plat[2][2]){
  this->puncture(1, plat[0][0], plat[1][0]);
  this->puncture(2, plat[0][1], plat[1][1]);
};



template<typename Type,typename PType>
void CubicIdeal<Type,PType> :: make_triangular(){

  Type bezout;
  XGCD(temp_mat[2][2], ci_temp, ci_temp2, coeff_matrix[2][0], coeff_matrix[2][1]);
  XGCD(temp_mat[2][2], bezout, ci_temp3, temp_mat[2][2], coeff_matrix[2][2]);
  mul(ci_temp, bezout, ci_temp);
  mul(ci_temp2, bezout, ci_temp2);

  // at this point, temp_mat[2][2] is the gcd of the bottom row
  // obtained by row operations rCol1 + sCol2 + tCol3

  mul(temp_mat[1][2], ci_temp, coeff_matrix[1][0]);
  mul(bezout, ci_temp2, coeff_matrix[1][1]);
  add(temp_mat[1][2], temp_mat[1][2], bezout);
  mul(bezout, ci_temp3, coeff_matrix[1][2]);
  add(temp_mat[1][2], temp_mat[1][2], bezout);

  mul(temp_mat[0][2], ci_temp, coeff_matrix[0][0]);
  mul(bezout, ci_temp2, coeff_matrix[0][1]);
  add(temp_mat[0][2], temp_mat[0][2], bezout);
  mul(bezout, ci_temp3, coeff_matrix[0][2]);
  add(temp_mat[0][2], temp_mat[0][2], bezout);

  // zero out the bottom row; affects the other entries
  coeff_matrix[1][0] -= (coeff_matrix[2][0]/temp_mat[2][2])*temp_mat[1][2];
  coeff_matrix[1][1] -= (coeff_matrix[2][1]/temp_mat[2][2])*temp_mat[1][2];
  coeff_matrix[1][2] -= (coeff_matrix[2][2]/temp_mat[2][2])*temp_mat[1][2];

  coeff_matrix[0][0] -= (coeff_matrix[2][0]/temp_mat[2][2])*temp_mat[0][2];
  coeff_matrix[0][1] -= (coeff_matrix[2][1]/temp_mat[2][2])*temp_mat[0][2];
  coeff_matrix[0][2] -= (coeff_matrix[2][2]/temp_mat[2][2])*temp_mat[0][2];

  XGCD(temp_mat[1][1], ci_temp, ci_temp2, coeff_matrix[1][0], coeff_matrix[1][1]);
  XGCD(temp_mat[1][1], bezout, ci_temp3, temp_mat[1][1], coeff_matrix[1][2]);
  mul(ci_temp, bezout, ci_temp);
  mul(ci_temp2, bezout, ci_temp2);


  mul(temp_mat[0][1], ci_temp, coeff_matrix[0][0]);
  mul(bezout, ci_temp2, coeff_matrix[0][1]);
  add(temp_mat[0][1], temp_mat[0][1], bezout);
  mul(bezout, ci_temp3, coeff_matrix[0][2]);
  add(temp_mat[0][1], temp_mat[0][1], bezout);


  coeff_matrix[0][0] -= (coeff_matrix[1][0]/temp_mat[1][1])*temp_mat[0][1];
  coeff_matrix[0][1] -= (coeff_matrix[1][1]/temp_mat[1][1])*temp_mat[0][1];
  coeff_matrix[0][2] -= (coeff_matrix[1][2]/temp_mat[1][1])*temp_mat[0][1];


  XGCD(temp_mat[0][0], ci_temp, ci_temp2, coeff_matrix[0][0], coeff_matrix[0][1]);
  XGCD(temp_mat[0][0], bezout, ci_temp3, temp_mat[0][0], coeff_matrix[0][2]);

  for (int i = 0; i<3; ++i){
    for (int j = 0; j<3; ++j){
      if (j >= i){
          coeff_matrix[i][j] = temp_mat[i][j];
      }else{
        ::clear(coeff_matrix[i][j]);
      }
    }
  }
  this->normalize();
};





template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: normalize(){

  this->ci_temp = this->denom;
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j<3; ++j){
      this->ci_temp = GCD(this->ci_temp, coeff_matrix[i][j]);
      if (IsOne(this->ci_temp) ){
        j =3; i = 3;
      }
    }
  }
    if (!IsOne(this->ci_temp) ){
      for (int i = 0; i < 3; ++i){
        for (int j = 0; j < 3; ++j){
          div(this->coeff_matrix[i][j],this->coeff_matrix[i][j],this->ci_temp );
        }
      }
    div(this->denom, this->denom, this->ci_temp);
    }
} // close normalize



template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: make_canonical(){

  this->normalize();

if ( !(IsZero(coeff_matrix[1][0]) && IsZero(coeff_matrix[2][0]) && IsZero(coeff_matrix[2][1]) ) )
{
  this->make_triangular();
}

    Type g, r, s;          // store gcd g and integer multipliers r,s which
                                 // satisfy r*m_32 + s*m_33 = g

    Type dummy, j,k,       //dummy is extraneous, j,k satisfy jr-ks = 1
         q,                // transitory value
         minorDet;         // determinant of minor matrix

    Type temp_mat[3][3]; //variable for the canonical matrix


    this->normalize();             // function which removes any common factors
                                  // in the denominator and numerators of entries

    // This if statement checks whether the properties of canonical basis are already
    //    satisfied

      if(
        ( this->coeff_matrix[2][1] == 0 ) &&
        ( 0 <= this->coeff_matrix[0][1] && this->coeff_matrix[0][1] < this->coeff_matrix[0][0]) &&
        ( 0 <= this->coeff_matrix[0][2] && this->coeff_matrix[0][2] < this->coeff_matrix[0][0]) &&
        ( 0 <= this->coeff_matrix[1][2] && this->coeff_matrix[1][2] < this->coeff_matrix[1][1])
      ) //close if clause
      {
        // do nothing
      }
      else{

      //computes the determinant of the lower-right 2x2 minor
      minorDet = (this->coeff_matrix[1][1] * this->coeff_matrix[2][2]
                  - this->coeff_matrix[1][2] * this->coeff_matrix[2][1]);

      //From the NTL Library, EEA
      XGCD(g, r, s ,this->coeff_matrix[2][1], this->coeff_matrix[2][2]);

      XGCD(dummy, j,k,r,-s);    // jr - ks = 1, dummy is always 1

      q = -(k*this->coeff_matrix[2][1] + j*this->coeff_matrix[2][2]);
      q = q/g;         // computes transitory value

      dummy = k+r*q;  //reusing dummy to save a couple ops


      temp_mat[0][1] =  (dummy) * this->coeff_matrix[0][1] + (j+q*s)*this->coeff_matrix[0][2];
      temp_mat[1][1] =  -minorDet/g;  // denoted  -epsilon/g in cubic book
      temp_mat[2][1] =  0;
      temp_mat[0][2] = r * this->coeff_matrix[0][1] +  s*this->coeff_matrix[0][2];
      temp_mat[1][2] = r * this->coeff_matrix[1][1] +  s*this->coeff_matrix[1][2];
      temp_mat[2][2] = g;

      // if -minorDet is negative, we need to flip the middle row to get e/g, where
      //  e is the absolute value of minorDet
      if (temp_mat[1][1] < 0){
          temp_mat[1][1] = -temp_mat[1][1];
          temp_mat[0][1] = - temp_mat[0][1];
      }
//std::cout << temp_mat[0][1]  << " " << temp_mat[0][2] << std::endl;
//std::cout << temp_mat[1][1]  << " " << temp_mat[1][2] << std::endl;
//std::cout << temp_mat[2][1]  << " " << temp_mat[2][2] << std::endl;
      //Ensures that temp_mat[1][1] (this is m_23 in the cubic book) satisfies
      // 0 <= m_23 < e/g )
      while ( (temp_mat[1][2] >= temp_mat[1][1]) || (temp_mat[1][2] < 0) ){
          if (temp_mat[1][2] < 0){
              temp_mat[0][2] += temp_mat[0][1];
              temp_mat[1][2] += temp_mat[1][1];
          }
          else if (temp_mat[1][2] >= temp_mat[1][1]){
              temp_mat[0][2] -= temp_mat[0][1];
              temp_mat[1][2] -= temp_mat[1][1];
          }
      }

      //ensures that m_12 is in the appropriate range
      while ( (temp_mat[0][1] >= this->coeff_matrix[0][0]) || (temp_mat[0][1] < 0) ){
          if (temp_mat[0][1] < 0){
              temp_mat[0][1] += this->coeff_matrix[0][0];
          }
          else if (temp_mat[0][1] >= this->coeff_matrix[0][0]){
              temp_mat[0][1] -= this->coeff_matrix[0][0];
          }
      }

      //ensures m_13 is in the appropriate range
      while ( (temp_mat[0][2] >= this->coeff_matrix[0][0]) || (temp_mat[0][2] < 0) ){
          if (temp_mat[0][2] < 0){
              temp_mat[0][2] += this->coeff_matrix[0][0];
          }
          else if (temp_mat[0][2] >= this->coeff_matrix[0][0]){
              temp_mat[0][2] -= this->coeff_matrix[0][0];
          }
      }

      for (int j = 1; j < 3 ; ++j){
        for (int i = 0; i <3; ++i){
          this->coeff_matrix[i][j] = temp_mat[i][j];
        }
      }
      this->normalize();
    }
} // close defn of make_canonical


template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: make_primitive(){
  ::clear(this->ci_temp);

  for (int j = 1; j < 3 ; ++j){
    for (int i = 0; i <3; ++i) {
      this->ci_temp = NTL::GCD(this->ci_temp, coeff_matrix[i][j]);
    }
  }

  if(!IsOne(this->ci_temp)){
    for (int j = 1; j < 3 ; ++j){
      for (int i = 0; i <3; ++i) {
        div(this->ci_temp,this->ci_temp, this->ci_temp);
      }
    }
  }
}


// By the end, adjacent_ideal will be the next ideal, while adj_min will contain the min adjacent to 1 of the last lattice
template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: adjacent_ideal(CubicIdeal<Type, PType> &B, CubicElement<Type, PType> &adj_min,char axis){

// Do I need to make sure the ideal is primitive? What about reduced?
  this->make_voronoi_basis(axis);

  this->divide_adjacent(B, adj_min);

  B.make_canonical();

}


// this function is intended to divide an ideal by the adjacent minima to 1,
// Assumes this ideal is already a Voronoi basis
template<typename Type, typename PType>
void CubicIdeal<Type, PType> :: divide_adjacent(CubicIdeal<Type, PType> &B, CubicElement<Type, PType> &adj_min){

    adj_min.set_order( this->my_order );
    adj_min.assign(this->get_coeff(0,1),this->get_coeff(1,1),this->get_coeff(2,1), this->get_coeff(0,0));

    // spare_ideal_element is set to be 1/theta_g
    this->spare_ideal_element.set_order(this->my_order);
    adj_min.inverse(this->spare_ideal_element);


    // set this->gen to be theta_h
    this->gen.set_order(this->my_order);
    this->gen.assign(this->get_coeff(0,2),this->get_coeff(1,2),this->get_coeff(2,2), this->get_coeff(0,0));

    std::cout << "DEBUG:"<<DEBUG << std::endl; 
    #ifdef DEBUG
    std::cout << "adjacent_ideal: The adjacent min is: "<< std::endl;
    std::cout << adj_min.get_u()  << " " << adj_min.get_x() << " " << adj_min.get_y() << std::endl;
    std::cout << "Below is the computed inverse: " << std::endl;
    std::cout << this->spare_ideal_element.get_u()  << " " << this->spare_ideal_element.get_x() << " " << this->spare_ideal_element.get_y() << std::endl;
    std::cout << "  "<< std::endl;
    CubicElement<Type, PType> checkValue = CubicElement<Type, PType>(this->my_order);
    mul(checkValue, adj_min, this->spare_ideal_element);
    std::cout << "inverse verification: Below should be '1 0 0'" << std::endl;
    std::cout << checkValue.get_u() << checkValue.get_x() << checkValue.get_y() << std::endl;
    #endif

    // multiplying with 1/theta_g yields theta_h/theta_g
    mul(this->gen, this->gen, this->spare_ideal_element);

    // Compute the LCM of the denominators
    this->ci_temp = GCD(this->gen.get_denom(), this->spare_ideal_element.get_denom());
    mul(this->ci_temp2,this->gen.get_denom(), this->spare_ideal_element.get_denom() );
    div(this->ci_temp, this->ci_temp2, this->ci_temp);


    // multiply our two elements by sigma = lcm (denom1, denom2)
    mul(this->gen, this->gen, this->ci_temp);
    mul(this->spare_ideal_element, this->spare_ideal_element, this->ci_temp);

    // Populates B to correspond to the lattice (1, 1/theta_g, theta_g/theta_h)* sigma
    B.assign( this->ci_temp, Type(0),Type(0),this->denom,
              this->spare_ideal_element.get_u(),this->spare_ideal_element.get_x(),this->spare_ideal_element.get_y(),this->denom,
              this->gen.get_u(),this->gen.get_x(),this->gen.get_y(),this->denom );


    #ifdef DEBUG
      std::cout << "Adj: The adjacent lattice below: "<< std::endl;
      std::cout << this->ci_temp  << " " << this->spare_ideal_element.get_u() << " " << this->gen.get_u() << std::endl;
      std::cout << Type(0)  << " " << this->spare_ideal_element.get_x() << " " << this->gen.get_x() << std::endl;
      std::cout << Type(0)  << " " << this->spare_ideal_element.get_y() << " " << this->gen.get_y() << std::endl;
    #endif


}


template<typename Type,typename PType>
bool CubicIdeal<Type, PType> :: is_equivalent(const CubicIdeal<Type, PType> & B){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: integral_norm( Type & val){
  // checks if the matrix is triangular
  if(( coeff_matrix[2][1] == Type(0)) && ( coeff_matrix[1][0] == Type(0)) && ( coeff_matrix[2][0] == Type(0)) ){
    mul(val, coeff_matrix[0][0], coeff_matrix[1][1]);
    mul(val, val, coeff_matrix[2][2]);
  }else{

    mul(this->ci_temp, coeff_matrix[1][1], coeff_matrix[2][2]);
    mul(this->ci_temp2, coeff_matrix[1][2],coeff_matrix[2][1]);
    sub(this->ci_temp, this->ci_temp, this->ci_temp2);
    mul(val, this->ci_temp, coeff_matrix[0][0]);


  //val.x = s23*y - x*s33;
    mul(this->ci_temp, coeff_matrix[1][2], coeff_matrix[2][0]);
    mul(this->ci_temp2, coeff_matrix[2][2],coeff_matrix[1][0]);
    sub(this->ci_temp, this->ci_temp, this->ci_temp2);
    mul(this->ci_temp, this->ci_temp, coeff_matrix[0][1]);
    add(val, val, this->ci_temp);

    //val.y = x*s32 - y*s22;
    mul(this->ci_temp, coeff_matrix[1][0], coeff_matrix[2][1]);
    mul(this->ci_temp2, coeff_matrix[1][1],coeff_matrix[2][0]);
    sub(this->ci_temp, this->ci_temp, this->ci_temp2);
    mul(this->ci_temp, this->ci_temp, coeff_matrix[0][2]);
    add(val, val, this->ci_temp);
  }
}
//s11*(s22*s33- s23*s32) + s12*(s12*s20- s22*s10) + s02*(s10*s21 - s11*s20)
template<typename Type,typename PType>
bool CubicIdeal<Type, PType> :: is_principal(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
bool CubicIdeal<Type, PType> :: is_prime(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
bool CubicIdeal<Type, PType> :: is_canonical(){
  std::cout <<"undefined" << std::endl;
  return false;
}

template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: reduce(CubicIdeal<Type, PType> & rIdeal,CubicElement<Type, PType> & relatingElement ){

  // relating element is to hold the element which you need to multiply the reduced ideal with the get back to the original
  // initialize by setting its order to this->ideal's, and assigning its value to 1.
  relatingElement.set_order(this->get_order());
  relatingElement.assign(Type(1));

  // copy the current element, rIdeal will hold the reduced ideal
  rIdeal.assign(*this);
  rIdeal.set_reduction_state(2);

  // This is going to hold an element which sits inside the normed body of 1
  CubicElement<Type, PType> interloper(this->get_order());

  while (rIdeal.is_reduced() != char(1)){
    // typical voronoi approach for reduction: Obtain a  Voronoi and divide by the adjacent minima
    // the behaviour of the make_voronoi algorithm is dependent on the ideal objects flag
    rIdeal.make_voronoi_basis();

    // because the is_reduced flag is not 1, the element in column 1 is an element
    // inside the norm body of 1. Set interloper to be this
    interloper.assign(rIdeal.get_coeff(0,1),rIdeal.get_coeff(1,1),rIdeal.get_coeff(2,1), rIdeal.get_coeff(0,0));



    rIdeal.divide_adjacent(rIdeal, interloper);
    mul(relatingElement, relatingElement,interloper);
    rIdeal.make_canonical();

    // for testing, print the real value of this element
    relatingElement.get_real_value(p_temp);
    //std::cout << " Ideal reduction: adj_min value before neg: "<< p_temp<< std::endl;
    if (p_temp < to<PType>(0)){
      relatingElement.negate(relatingElement);
    }
      //relatingElement.get_real_value(p_temp);
    //std::cout << " Ideal reduction: adj_min value: "<< p_temp<< std::endl;
  }
}


template<typename Type,typename PType>
void CubicIdeal<Type, PType> :: make_voronoi_basis(char axis){
  #ifdef DEBUGVORONOI
  std::cout << "Computing Voronoi basis..." << std::endl;
  #endif
  if (axis == 'X'){
    this->get_order()->get_voronoi()->make_voronoi_basis(*this);
  }
  else if (axis == 'Z'){
    this->get_order()->roots_swap_position(0,2);
    this->get_order()->get_voronoi()->make_voronoi_basis(*this);
    this->get_order()->roots_swap_position(0,2);

  }else if (axis == 'Y'){
    this->get_order()->roots_swap_position(0,1);
    this->get_order()->get_voronoi()->make_voronoi_basis(*this);
    this->get_order()->roots_swap_position(0,1);
  }




}


/// Friend definitions
// This function might modify the ideal representation!
template<typename T,typename PT>
bool is_equal(CubicIdeal<T,PT> & A, CubicIdeal <T,PT> & B){
  #ifdef DEBUG
  std::cout << "... Comparing Lattices ... "<< std::endl;
  #endif
  /*LatticeBasis B1,B2;
  for (int i = 0; i<3; ++i){
    for (int j = 0; j<2; ++j){
      B1.coefficientMatrix[i][j] = L1.coefficientMatrix[i][j];
      B2.coefficientMatrix[i][j] = L2.coefficientMatrix[i][j];
    }
  }
  B1.mainDenominator = L1.mainDenominator;
  B.coeff_matrix[0][0] = L2.mainDenominator;


    /*potentially add an if clause that checks whether a basis is already
    canonical. This would require a latticeBasis object to have a bool that
    tells us */

    A.make_canonical();
    B.make_canonical();

    // store the values e_1 and e_2 in static variables
    mul(A.ci_temp,  A.coeff_matrix[1][1], A.coeff_matrix[2][2]);
    mul(A.ci_temp3, A.coeff_matrix[1][2], A.coeff_matrix[2][1]);
    sub(A.ci_temp,  A.ci_temp, A.ci_temp3);
    abs(A.ci_temp,  A.ci_temp);
    // A.ci_temp = abs(A.coeff_matrix[1][1] * A.coeff_matrix[2][2] - A.coeff_matrix[1][2] * A.coeff_matrix[2][1]);
    mul(A.ci_temp2, B.coeff_matrix[1][1], B.coeff_matrix[2][2]);
    mul(A.ci_temp3, B.coeff_matrix[1][2], B.coeff_matrix[2][1]);
    sub(A.ci_temp2, A.ci_temp2, A.ci_temp3);
    abs(A.ci_temp2, A.ci_temp2);
    //A.ci_temp2 = abs(B.coeff_matrix[1][1] * B.coeff_matrix[2][2] - B.coeff_matrix[1][2] * B.coeff_matrix[2][1]);


    if (   ( A.coeff_matrix[0][0] != B.coeff_matrix[0][0]   )             // matching sigma values
        || ( A.ci_temp != A.ci_temp2 )                            // e values match
        || (A.coeff_matrix[2][2] != B.coeff_matrix[2][2]) ){              //g values match

        #ifdef DEBUG
        std::cout << "first lattice-equality check fails" << std::endl;
        #endif
        return false;
    }//close if clause

    else if (
         (  ( (A.coeff_matrix[1][2] - B.coeff_matrix[1][2]) % A.coeff_matrix[1][1] ) != 0  )
//      || ( (A.coeff_matrix[0][1]%A.coeff_matrix[0][0]) != (B.coeff_matrix[0][1] % A.coeff_matrix[0][0]) )
      || (  ( (A.coeff_matrix[0][1] - B.coeff_matrix[0][1]) %A.coeff_matrix[0][0]) != 0 )
      || (   ( (A.coeff_matrix[0][2] - B.coeff_matrix[0][2])%A.coeff_matrix[0][0])
          != ( ( (A.coeff_matrix[1][2] - B.coeff_matrix[1][2])*A.coeff_matrix[0][1] / A.coeff_matrix[1][1] ) % A.coeff_matrix[0][0])   )
      ) //close else if clause
    {
        #ifdef DEBUG
        std::cout << "second latice-equality check fails" << std::endl;
        #endif

        return false;
    }
    else if (A.get_denom() != B.get_denom()){
      return false;       // added condition after adapting since now we might have fractional ideals
    }
    else {
      return true;
    }
} //close function compareLattice3


template <typename Type, typename PType>
void mul(CubicIdeal<Type,PType> & A, const CubicIdeal<Type,PType> & B, const CubicIdeal<Type,PType> & C){
  if (!(A.get_order()->is_equal( (*B.get_order()) ) ) || (!(A.get_order()->is_equal( (*C.get_order()) ) ) )){
    std::cout << "Order mismatch in ideal mul" << std::endl;
  }
  (A.get_order()->m_method)->multiply(A,B,C);
}














#endif
