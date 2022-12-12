#include <ANTL/Quadratic/Reduce/ReducePlainReal_Opt.hpp>

// ReducePlainReal<ZZ>::Reduce
// Task: Reduces the ideal

template <> void ReducePlainRealOpt<ZZ>::reduce(QuadraticIdealBase<ZZ> &A) {

  static ZZ q, r, temp, a2, nb, na, s;
  static ZZ a, b, c, rootD;
  static ZZ temp_num_q;

  a = A.get_a();
  b = A.get_b();
  c = A.get_c();
  rootD = FloorRootDelta;

  A.set_num_q(ZZ(0));

  // na = (D - b^2) / 4a
  c = -c;

  a2 = abs(a) << 1;
  temp = rootD - a2;
  if (temp < 0)
    ++temp;

  while ((abs(temp) >= b) || (b > rootD) || (a < 0)) {
    s = (a > 0) ? 0 : 1;

    // (rootD+b) = (2a)q + r
    a2 = a << 1;
    temp = rootD+b;
    if (to_long(s))  ++temp;
    DivRem(q,r,temp,a2);

    temp_num_q = A.get_num_q();
    A.set_qlist_i(to_long(temp_num_q), q);
    A.set_num_q(temp_num_q + ZZ(1));

    // nb = rootD + s - r;
    nb = rootD - r;
    if (to_long(s))  ++nb;

    // na = c - q * (nb - b)/2
    na = c - q*((nb - b) >> 1);

    b = nb;
    c = a;
    a = na;

    a2 = abs(a) << 1;
    temp = rootD - a2;
    if (temp < 0)
       ++temp;
  }

  c = -c;
  //  assign_abc (qie, a, b, c);
  A.assign(a, b, c);
}

