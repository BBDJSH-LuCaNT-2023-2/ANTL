#ifndef ANTL_CUBIC_ELEMENTNF_CPP
#define ANTL_CUBIC_ELEMENTNF_CPP

#include "../../include/ANTL/Cubic/CubicElementNF.hpp"


template class CubicElementNF<long, double>;
// static member definitions

template<typename Type, typename PType>
Type CubicElementNF<Type, PType>::temp2;


template <typename Type, typename PType>
CubicElementNF<Type, PType> :: CubicElementNF(const CubicOrderNF<Type, PType> * cnfo, const Type _coefficients[3], const Type & denom)
  : CubicElement<Type,PType>::CubicElement(cnfo, _coefficients, denom){

    reduce();
  }

/*
template <typename Type, typename PType>
void CubicElementNF<Type, PType> :: norm(ANTL::QQ<Type> & newVal){

  boost::math::tools::polynomial<Type> order_poly = (*this->my_order).get_IBCF();

  Type a,b,c,d,u,x,y;
  u = this->coefficients[0]; x = this->coefficients[1]; y = this->coefficients[2];
  a = order_poly[3]; b = order_poly[2];
  c = order_poly[1]; d = order_poly[0];
  newVal.setNumerator( (power(u,3) - power(u,2)*x*b - 2*c*power(u,2)*y + a*c*u*power(x,2)
      +(b*c+ 3*a*d)*u*x*y + (b*d + power(c,2))*u*power(y,2)
      - power(a,2)*d*power(x,3) - 2*a*b*d*power(x,2)*y
      - d*(power(b,2)+a*c)*x*power(y,2) +d*(a*d -b*c)*power(y,3) ));


  //newVal.setNumerator( (power(u,3) - power(u,2)*x*b - 2*c*power(u,2)*y + a*c*u*power(x,2)
  //    +(b*c+ 3*a*d)*u*x*y + (b*d + power(c,2))*u*power(y,2)
  //    - power(a,2)*d*power(x,3) - 2*a*b*d*power(x,2)*y
  //    - d*(power(b,2)+a*c)*x*power(y,2) +d*(a*d -b*c)*power(y,3) ));

  newVal.setDenominator(this->denom^3 ) ;
}

*/


//
template <typename Type, typename PType>
void CubicElementNF<Type, PType> :: trace(ANTL::QQ<Type> & newVal){
  boost::math::tools::polynomial<Type> order_poly = (*this->my_order).get_IBCF();
  newVal.setNumerator( 3*this->u
  - order_poly[2] *this->x
  - 2 * order_poly[1] * this->y);

  newVal.setDenominator(this->denom);
}

/*
template <typename Type, typename PType>
void CubicElementNF<Type, PType> :: inverse(CubicElement<Type, PType> & val){
  std::cout << "TBI "<< std::endl;
  // Using the multiplication table, we need to find r,s,t so that
  // (r + s*p1 + t*p2)*(u + x*p1 + y*p2) = 1
  // Multiplying out , we get a system of equations
  // ur +       (-ady)s +  (-adx -bdy)t = 1
  // xr + (u - bx -cy)s +   (-cx - dy)t = 0
  // yr +         (ax)s +       (u-cy)t = 0
  //That can be represented by a matrix S

  //Solve the matrix S with Cramers Rule

  Type u,x, y;
  Type a,b,c,d;
  Type s11, s12, s13, s21, s22, s23, s31, s32,s33;
  u = this->u;
  x = this->x;
  y = this->y;
  a = (this->my_order)->get_IBCF()[3];
  b = (this->my_order)->get_IBCF()[2];
  c = (this->my_order)->get_IBCF()[1];
  d = (this->my_order)->get_IBCF()[0];

  s11 = u;
  s12 = - a*d*y;
  s13 = -a*d*x - b*d*y;
  s21 = x;
  s22 = u - b*x-c*y;
  s23 = -c*x -d*y;
  s31 = y;
  s32 = a*x;
  s33 = u - c*y;

  Type Disc = s11*s22*s33 + s12*s23*s31 + s13*s21*s32 - s23*s32*s11 - s33*s21*s12 - s13*s22*s31;

  //val.u = s22*s33 - s23*s32;
  //val.x = s23*s31 - s21*s33;
  //val.y = s21*s32 - s31*s22;

}*/



template <typename Type, typename PType>
bool CubicElementNF<Type, PType> :: is_zero() const {
  return (this->u ==Type(0)) && (this->x ==Type(0)) && (this->y == Type(0));
}

template <typename Type, typename PType>
bool CubicElementNF<Type, PType> :: is_equal(const CubicElementNF <Type, PType> & B){
  return ( ( (this->my_order)->is_equal(*B.my_order) ) && (this->u == B.u) && (this->x == B.x) && (this->y == B.y)&& (this->denom == B.denom));
}


template <typename Type, typename PType>
void CubicElementNF<Type, PType> :: assign(const CubicElementNF<Type,PType> & C){
  this->u = C.u;
  this->x = C.x;
  this->y = C.y;
  this->denom = C.denom;

}



/****************************************************************************
*************************** friend functions ********************************/


template <typename T, typename PT>
void assign (CubicElement <T,PT> & A, const CubicElement <T,PT> & B){
  if (!(A.my_order.is_equal(B.my_order)) ){
    std::cout << "Orders are not equal, throw exception" << std::endl;
  }
  else{
    A.u = B.u;
    A.x = B.x;
    A.y = B.y;
    A.denom = B.denom;
  }
}
/****** Operations *******/


template <typename T, typename PT>
void add (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT> & C){


  //  (B.u * C.denom + C.u * B.denom) + rho1*(B.x * C.denom + C.x * B.denom ) + rho2*(B.y * C.denom + C.y * B.denom )
  //    -------------------------------------------------------------------------------------------------------------------
  //                                                    B.denom * C.denom
  //

  if (!( B.my_order->is_equal( *(C.my_order) ) ) ){
    std::cout << "element orders are not the same. Throw exception" << std::endl;
  }

  if (B.is_zero()){
    A.assign(C);
  }
  else if (C.is_zero()){
    A.assign(B);
  }
  else{
    mul(A.u, B.u, C.denom); // u1*d2
    mul(CubicElement<T,PT>::temp, C.u, B.denom);
    add(A.u, A.u, CubicElement<T,PT>::temp);

    mul(A.x, B.x, C.denom); // u1*d2
    mul(CubicElement<T,PT>::temp, C.x, B.denom);
    add(A.x, A.x, CubicElement<T,PT>::temp);

    mul(A.y, B.y, C.denom); // u1*d2
    mul(CubicElement<T,PT>::temp, C.y, B.denom);
    add(A.y, A.y, CubicElement<T,PT>::temp);

    A.reduce();
  }
}

/* Add a constant integer to the Cubic Number A */
template <typename T, typename PT>
void add (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const T & alpha){

  if (B.is_zero()){
    A.u = alpha;
    ::clear(A.x);
    ::clear(A.y);
    ::set(A.denom);
  }
  else if (IsZero(alpha)){
    A.assign(B);
  }
  else{
    mul(CubicElement<T,PT>::temp, alpha, B.denom);
    add(A.u, B.u, CubicElement<T,PT>::temp);

    A.x = B.x;
    A.y = B.y;
    A.denom = B.denom;
    A.reduce();
  }
}

template <typename T, typename PT>
void add (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const QQ<T> & alpha){
  if (!A.my_order.is_equal(B.my_order)){
    std::cout << "order mismatch, throw exception" << std::endl;
  }

  if ( B.is_zero() ){
    A.assign(alpha);
  }
  else if (IsZero(alpha)){
    A.assign(B);
  }
  else{
    mul(CubicElement<T,PT>::newU, B.u, alpha.getDenominator());
    mul(CubicElement<T,PT>::temp, alpha.getNumerator(), B.denom);
    add(CubicElement<T,PT>::newU, CubicElement<T,PT>::newA, CubicElement<T,PT>::temp);

    mul(CubicElement<T,PT>::newX, B.x, alpha.getDenominator());
    mul(CubicElement<T,PT>::newY, B.y, alpha.getDenominator());
    mul(CubicElement<T,PT>::temp, B.denom, alpha.getDenominator());

    A.u = CubicElement<T,PT>::newU;
    A.x = CubicElement<T,PT>::newX;
    A.y = CubicElement<T,PT>::newY;
    A.denom = CubicElement<T,PT>::temp;

    A.reduce();

  }
}

// Multiplication
// see CFG pg 40
template <typename T, typename PT>
void mul (CubicElementNF <T,PT> & A, const CubicElementNF <T,PT> & B, const CubicElementNF <T,PT> & C){


  // COmpute Y
  mul(CubicElement<T,PT>::newU, C.x,B.y);                                                           // x_2* y_1
  mul(CubicElement<T,PT>::temp, C.y, B.x);                                                          // x_1* y_2
  add(CubicElement<T,PT>::newU, CubicElement<T,PT>::newU, CubicElementNF<T, PT>::temp);             // x_2* y_1 + x_1* y_2

  mul(CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[0] ,(A.my_order)->defining_IBCF[3] );   // ad
  mul(CubicElement<T,PT>::newU, -CubicElement<T,PT>::temp, CubicElement<T,PT>::newU);               // -ad * (x2*y_1 + x_1*y_2)
  mul(CubicElement<T,PT>::temp, C.u, B.u);                                                          // u1*u2
  add(CubicElement<T,PT>::newU, CubicElement<T,PT>::newU, CubicElement<T,PT>::temp);                // u1*u2 - ad*(x2*y_1 + x_1*y_2)

  mul(CubicElement<T,PT>::temp, C.y, B.y);                                                          // y1*y2
  mul(CubicElement<T,PT>::temp, CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[3]);          // d*y1*y2
  mul(CubicElement<T,PT>::temp, CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[1]);          // b*d*y1*y2

  sub(CubicElement<T,PT>::newU, CubicElement<T,PT>::newU, CubicElement<T,PT>::temp);                // u1*u2 - ad*(x2*y_1 + x_1*y_2) - b*d*y1*y2
  //mul(temp, temp, );

  // Compute X
  mul(CubicElement<T,PT>::newX, B.x,C.y);                                                           // x1*y2
  mul(CubicElement<T,PT>::temp, C.x,B.y);                                                           // x2*y1
  add(CubicElement<T,PT>::newX, CubicElement<T,PT>::newX,CubicElement<T,PT>::temp);                 // x1*y2 + x2*y1
  mul(CubicElement<T,PT>::newX, CubicElement<T,PT>::newX, -(A.my_order)->defining_IBCF[2]);         // -c*(x1*y2 + x2*y1)

  mul(CubicElement<T,PT>::temp, B.u, C.x);                                                          // u1*x2
  add(CubicElement<T,PT>::newX, CubicElement<T,PT>::newX, CubicElement<T,PT>::temp);                // u1*x2 -c*(x1*y2 + x2*y1)

  mul(CubicElement<T,PT>::temp, B.x, C.u);                                                          // u2*x1
  add(CubicElement<T,PT>::newX, CubicElement<T,PT>::newU, CubicElement<T,PT>::temp);                // u1*x2 + u2*x1 -c*(x1*y2 + x2*y1)

  mul(CubicElement<T,PT>::temp, C.x, B.x);                                                          // x1*x2
  mul(CubicElement<T,PT>::temp, CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[1]);          // b*x1*x2
  sub(CubicElement<T,PT>::newX, CubicElement<T,PT>::newX, CubicElement<T,PT>::temp);                // u1*x2 + u2*x1 -b*x1*x2 - c*(x1*y2 + x2*y1)

  mul(CubicElement<T,PT>::temp, C.y, B.y);                                                          // y1*y2
  mul(CubicElement<T,PT>::temp, CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[3]);          // d*y1*y2
  sub(CubicElement<T,PT>::newX, CubicElement<T,PT>::newX, CubicElement<T,PT>::temp);                // u1*x2 + u2*x1 -b*x1*x2 - c*(x1*y2 + x2*y1) - d*y1*y2

  // Compute X
  mul(CubicElement<T,PT>::newY, B.u,C.y);                                                           // u1*y2
  mul(CubicElement<T,PT>::temp, C.u,B.y);                                                           // u2*y1
  add(CubicElement<T,PT>::newY, CubicElement<T,PT>::newY,CubicElement<T,PT>::temp);                 // u1*y2 + u2*y1

  mul(CubicElement<T,PT>::temp, C.x, B.x);                                                          // x1*x2
  mul(CubicElement<T,PT>::temp, CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[0]);          // a*x1*x2
  add(CubicElement<T,PT>::newY, CubicElement<T,PT>::newY, CubicElement<T,PT>::temp);                // u1*y2 + u2*y1 + a*x1*x2

  mul(CubicElement<T,PT>::temp, C.y, B.y);                                                          // y1*y2
  mul(CubicElement<T,PT>::temp, CubicElement<T,PT>::temp, (A.my_order)->defining_IBCF[2]);          // c*y1*y2
  sub(CubicElement<T,PT>::newY, CubicElement<T,PT>::newY, CubicElement<T,PT>::temp);                // u1*y2 + u2*y1 + a*x1*x2 - c*y1*y2


  mul(A.denom, C.denom, B.denom);

  A.u = CubicElement<T,PT>::newU;
  A.x = CubicElement<T,PT>::newX;
  A.y = CubicElement<T,PT>::newY;

  A.reduce();

}

template <typename T, typename PT>
void sub (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const CubicElement <T,PT> & C){
  if (!( (B.get_order())->is_equal( *(C.get_order())))){
    std::cout << "element orders are not the same. Throw exception" << std::endl;

  }
  else{
    mul(A.u, B.u, C.denom); // u1*d2
    mul(CubicElement<T,PT>::temp, C.u, B.denom);
    sub(A.u, A.u, CubicElement<T,PT>::temp);

    mul(A.x, B.x, C.denom); // u1*d2
    mul(CubicElement<T,PT>::temp, C.x, B.denom);
    sub(A.x, A.x, CubicElement<T,PT>::temp);

    mul(A.y, B.y, C.denom); // u1*d2
    mul(CubicElement<T,PT>::temp, C.y, B.denom);
    sub(A.y, A.y, CubicElement<T,PT>::temp);

    A.reduce();
  }
}

/* Sub a constant integer to the Cubic Number A */
template <typename T, typename PT>
void sub (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const T & alpha){

  if (B.is_zero()){
    A.u = -alpha;
    ::clear(A.x);
    ::clear(A.y);
    ::set(A.denom);
  }
  else if (IsZero(alpha)){
    A.assign(B);
  }
  else{
    mul(CubicElement<T,PT>::temp, alpha, B.denom);
    sub(A.u, B.u, CubicElement<T,PT>::temp);

    A.x = B.x;
    A.y = B.y;
    A.denom = B.denom;
    A.reduce();
  }
}

template <typename T, typename PT>
void sub (CubicElement <T,PT> & A, const CubicElement <T,PT> & B, const QQ<T> & alpha){
  if (!A.my_order.is_equal(B.my_order)){
    std::cout << "order mismatch, throw exception" << std::endl;
  }

  if ( B.is_zero() ){
    A.assign(-alpha);
  }
  else if (IsZero(alpha)){
    A.assign(B);
  }
  else{
    mul(CubicElement<T,PT>::newU, B.u, alpha.getDenominator());
    mul(CubicElement<T,PT>::temp, alpha.getNumerator(), B.denom);
    add(CubicElement<T,PT>::newU, CubicElement<T,PT>::newA, CubicElement<T,PT>::temp);

    mul(CubicElement<T,PT>::newX, B.x, alpha.getDenominator());
    mul(CubicElement<T,PT>::newY, B.y, alpha.getDenominator());
    mul(CubicElement<T,PT>::temp, B.denom, alpha.getDenominator());

    A.u = CubicElement<T,PT>::newU;
    A.x = CubicElement<T,PT>::newX;
    A.y = CubicElement<T,PT>::newY;
    A.denom = CubicElement<T,PT>::temp;

    A.reduce();

  }
}






template <typename Type, typename PType>
void divide (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const CubicElement <Type,PType> & C){

}

template <typename Type, typename PType>
void power (CubicElement <Type,PType> & A, const CubicElement <Type,PType> & B, const NTL::ZZ & p){

}
template <typename Type, typename PType>
void CubicElementNF<Type, PType> :: reduce(){
  std::cout << "implement reduction" << std::endl;
}



#endif
