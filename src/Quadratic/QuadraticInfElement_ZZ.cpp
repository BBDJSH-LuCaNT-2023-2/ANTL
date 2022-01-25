#include <ANTL/Quadratic/QuadraticInfElement.hpp>

using namespace ANTL;

template <>
void QuadraticInfElement<ZZ>::baby_step() {
  // A variable Xi refers to X_{i+i} E.g. R2 refers to R_{i+2}
  ZZ R, Q, P, B;
  ZZ a, b, c;
  ZZ q, r;

  RR RelativeGenerator;

  a = qib.get_a();
  b = qib.get_b();
  c = qib.get_c();

  DivRem(q, r, b + FloorRootDelta, 2*a);

  R = -a;
  P = FloorRootDelta - r;
  Q = q*((b - P)/2) - c;

  qib.assign(Q, P, R);
  qib.normalize();
  Distance += log(inv(abs((to_RR(P) - sqrt(to_RR(Delta))) / to_RR(2*Q))));
}

template <>
void QuadraticInfElement<ZZ>::giant_step(QuadraticInfElement<ZZ> & quad_ib) {

  mul(qib, qib, quad_ib.get_qib());

  Distance += quad_ib.get_distance();
  Distance += log(qib.get_QO()->get_mul_nucomp()->get_RelativeGenerator()->conv_RR());

}

//Debug tools
//   if(debug) {
//     std::cout << "a is " << a << std::endl;
//     std::cout << "b is " << b << std::endl;
//     std::cout << "c is " << c << std::endl;
//     std::cout << "q is " << q << std::endl;
//     std::cout << "r is " << r << std::endl;
//     std::cout << "Q is " << Q << std::endl;
//     std::cout << "P is " << P << std::endl;
//     std::cout << "R is " << R << std::endl;
//     std::cout << "RG is " << RelativeGenerator << std::endl;
//     std::cout << "distance is " << log(RelativeGenerator) << std::endl;
//     std::cout << "Distance is " << Distance << std::endl;
//   }
