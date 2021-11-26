#include <ANTL/Quadratic/QuadraticInfElement.hpp>

using namespace ANTL;

template <>
void QuadraticInfElement<ZZ>::baby_step() {
  // A variable Xi refers to X_{i+i} E.g. R2 refers to R_{i+2}
  bool debug = false;
  ZZ R, Q, P, B;
  ZZ a, b, c;
  ZZ q, r;

  RR RelativeGenerator;

  a = qib.get_a();
  b = qib.get_b();
  c = qib.get_c();

  DivRem(q, r, b + FloorRootDelta, 2*a);

  if(debug) {
    std::cout << "a is " << a << std::endl;
    std::cout << "b is " << b << std::endl;
    std::cout << "c is " << c << std::endl;
    std::cout << "q is " << q << std::endl;
    std::cout << "r is " << r << std::endl;
  }

  R = -a;
  P = FloorRootDelta - r;
  Q = q*((b - P)/2) - c;

  RelativeGenerator = inv(abs((to_RR(P) - sqrt(to_RR(Delta))) / to_RR(2*Q)));

  if(debug) {
    std::cout << "Q is " << Q << std::endl;
    std::cout << "P is " << P << std::endl;
    std::cout << "R is " << R << std::endl;
    std::cout << "RG is " << RelativeGenerator << std::endl;
    std::cout << "distance is " << log(RelativeGenerator) << std::endl;
  }


  qib.assign(Q, P, R);
  distance = log(RelativeGenerator);
}
