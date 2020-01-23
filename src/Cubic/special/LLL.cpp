#ifndef CUBIC_LLL_H
#define CUBIC_LLL_H
#include "../../include/ANTL/Cubic/special/LLL.hpp"

// I think the input should in unit rank 2 be a 2 by 2 matrix with real entries
const RR DOUBLE_TOL(1e-25);




template<PType> void LLL_reduction(PType logbasis[2][2]){
  LLL_reduction(logbasis[0][0], logbasis[1][0], logbasis[0][1], logbasis[1][1]);
}

template<PType> void LLL_reduction(PType lambda00, PType lambda10, PType lambda01, PType lambda11){
  int k = 1;
  int kmax = 0;
  PType mu10, LLLtemp, B1, q;
  PType bstar0[2] = {lambda00, lambda10};
  PType bstar1[2];

  PType B0 = bstar0[0]*bstar0[0] + bstar0[1]* bstar0[1];
  Ptype H[2][2] = {{PType(1), PType(0)},{PType(0), PType(1)}};


  if (k > kmax){
    // k = 1
    kmax = k;
    bstar1[0] = lambda01;
    bstar1[1] = lambda11;
    mul(mu10, bstar0[0],bstar1[0]);
    mul(LLLtemp, bstar0[1],bstar1[1]);
    add(mu10, LLLtemp, mu10);
    div(mu10,mu10, B0);
    B1 = bstar1[0]*bstar1[0] + bstar1[1]*bstar1[1];
    if (abs(B1) <= DOUBLE_TOL){
      std::cout << "Error, not a basis" << std::endl;
    }
  }
  bool flag;

  do{
    flag = true;
  // sub algorithm red
  if (abs(mu10) > 0.5){
    add(q, mu10, 0.5);
    floor(q, q);

    mul(LLLtemp, bstar0[0], q);
    sub(bstar1[0], bstar1[0], LLLtemp);
    mul(LLLtemp, bstar0[1], q);
    sub(bstar1[1], bstar1[1], LLLtemp);

    sub(H[0][1], H[0][1] - q);
    sub(mu10, mu10, q);
  }

  mul(LLLtemp, mu10, mu10);
  sub(LLLtemp, LLLtemp, 0.75);
  LLLtemp = -LLLtemp;
  mul(LLLtemp, B0);
  if(B1 < LLLtemp ){
    std::swap(bstar0[0], bstar1[0]);
    std::swap(bstar0[1], bstar1[1]);
    std::swap(H[0][0], H[0][1]);
    std::swap(H[1][0], H[1][1]);

    PType mu,B, b0,b1;
    mu = mu10;
    B = B1- (mu^2)*B0;
    mu10 = mu*B0/B;

    b0 = bstar0[0];
    b1 = bstar0[1];

    bstar0[0]=bstar1[0] + mu*b0;
    bstar0[1]=bstar1[1] + mu*b1;

    bstar1[0] = -mu10*bstar1[0] +(B1/B)*b0;
    bstar1[1] = -mu10*bstar1[1] +(B1/B)*b1;

    B1 = B0*B1/B;
    B0 = B;

    flag = false;
  }
  while(flag);

}










#endif // guard
