#ifndef ANTL_IDEAL_MULTIPLICATION_STRATEGY_CPP
#define ANTL_IDEAL_MULTIPLICATION_STRATEGY_CPP

#include "../../../include/ANTL/Cubic/Multiplication/IdealMultiplicationStrategy.hpp"

template<typename Type, typename PType>
void IdealMultiplicationStrategy<Type, PType> :: column_mul(Type & out1, Type & out2, Type & out3, const Type & u1, const Type & x1, const Type  & y1,
const Type & u2, const Type & x2, const Type & y2){

    //error check
    if (my_order == NULL){
      std::cout << "my_order is NULL " << std::endl;
      throw std::exception();
    }
    // Compute U
    mul(out1, x2,y1);                                                           // x_2* y_1
    mul(multemp, x1, y2);                                                          // x_1* y_2
    add(out1, out1, multemp);             // x_2* y_1 + x_1* y_2

    mul(multemp, (this->my_order)->get_coeff(0) ,(this->my_order)->get_coeff(3) );   // ad
    mul(out1, -multemp, out1);               // -ad * (x2*y_1 + x_1*y_2)
    mul(multemp, u2, u1);                                                          // u1*u2
    add(out1, out1, multemp);                // u1*u2 - ad*(x2*y_1 + x_1*y_2)

    mul(multemp, y2, y1);                                                          // y1*y2
    mul(multemp, multemp, (this->my_order)->get_coeff(0));          // d*y1*y2
    mul(multemp, multemp, (this->my_order)->get_coeff(2));          // b*d*y1*y2

    sub(out1, out1, multemp);                // u1*u2 - ad*(x2*y_1 + x_1*y_2) - b*d*y1*y2
    //mul(temp, temp, );

    // Compute X
    mul(out2, x1,y2);                                                           // x1*y2
    mul(multemp, x2,y1);                                                           // x2*y1
    add(out2, out2,multemp);                 // x1*y2 + x2*y1
    mul(out2, out2, -(this->my_order)->get_coeff(1));         // -c*(x1*y2 + x2*y1)

    mul(multemp, u1, x2);                                                          // u1*x2
    add(out2, out2, multemp);                // u1*x2 -c*(x1*y2 + x2*y1)

    mul(multemp, x1, u2);                                                          // u2*x1
    add(out2, out2, multemp);                // u1*x2 + u2*x1 -c*(x1*y2 + x2*y1)

    mul(multemp, x2, x1);                                                          // x1*x2
    mul(multemp, multemp, (this->my_order)->get_coeff(2));          // b*x1*x2
    sub(out2, out2, multemp);                // u1*x2 + u2*x1 -b*x1*x2 - c*(x1*y2 + x2*y1)

    mul(multemp, y2, y1);                                                          // y1*y2
    mul(multemp, multemp, (this->my_order)->get_coeff(0));          // d*y1*y2
    sub(out2, out2, multemp);                // u1*x2 + u2*x1 -b*x1*x2 - c*(x1*y2 + x2*y1) - d*y1*y2

    // Compute Y
    mul(out3, u1,y2);                                                           // u1*y2
    mul(multemp, u2,y1);                                                           // u2*y1
    add(out3, out3,multemp);                 // u1*y2 + u2*y1

    mul(multemp, x2, x1);                                                          // x1*x2
    mul(multemp, multemp, (this->my_order)->get_coeff(3));          // a*x1*x2
    add(out3, out3, multemp);                // u1*y2 + u2*y1 + a*x1*x2

    mul(multemp, y2, y1);                                                          // y1*y2
    mul(multemp, multemp, (this->my_order)->get_coeff(1));          // c*y1*y2
    sub(out3, out3, multemp);                // u1*y2 + u2*y1 + a*x1*x2 - c*y1*y2


    //mul(A.denom, C.denom, B.denom);

    //A.u = CubicElement<T,PT>::newU;
    //A.x = CubicElement<T,PT>::newX;
    //A.y = CubicElement<T,PT>::newY;

    //A.normalize();

}








#endif
