using namespace NTL;
using namespace ANTL;

template <class T> void QuadraticNumber<T>::normalize () {
  // remove common factors
  T g = GCD(GCD(a,b),d);
  if (!::IsOne(g)) {
      ::div(a,a,g);
      ::div(b,b,g);
      ::div(d,d,g);
  }

  // normalize leading coefficients
  MakeMonic(d);
  if (deg(a) > (deg(b) + QO->getGenus()))
    MakeMonic(a);
  else
    MakeMonic(b);
}
