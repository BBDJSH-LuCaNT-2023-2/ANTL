using namespace ANTL;

namespace ANTL {

template <class T, class S>
QuadraticInfElement<T, S>::QuadraticInfElement(QuadraticOrder<T> &quad_o)
    : qib{quad_o} {
  qib.assign_one();
  Delta = quad_o.get_discriminant();
  FloorRootDelta = FloorToZZ(sqrt(to_RR(Delta)));
  Distance = 0;
}

template <class T, class S> QuadraticInfElement<T, S>::~QuadraticInfElement() {}

template <class T, class S>
QuadraticIdealBase<T> QuadraticInfElement<T, S>::get_qib() const {
  return qib;
}

template <class T, class S> S QuadraticInfElement<T, S>::get_distance() const {
  return Distance;
}

// Infrastructure methods
template <class T, class S> void QuadraticInfElement<T, S>::baby_step() {
  T a = qib.get_a(), b = qib.get_b(), c = qib.get_c(), q, r, R, Q, P;
  S relative_generator, relative_distance;

  // Computation of the QuadraticIdealBase via a single step in the continued fraction expansion
  DivRem(q, r, b + FloorRootDelta, 2 * a);

  R = -a;
  P = FloorRootDelta - r;
  Q = q * ((b - P) / 2) - c;

  qib.assign(Q, P, R);
  qib.normalize();

  // Computation of the relative distance
  relative_generator = abs(to<S>(2 * Q) / (to<S>(P) - sqrt(to<S>(Delta))));
  relative_distance = log(relative_generator);

  // Updating the distance
  update_distance_add(Distance, Distance, relative_distance);
}

template <class T, class S>
void QuadraticInfElement<T, S>::giant_step(
    const QuadraticInfElement<T, S> &qif_1) {
  S relative_distance_1, relative_distance_2;

  mul(qib, qib, qif_1.get_qib());

  relative_distance_1 = qif_1.get_distance();
  relative_distance_2 = qib.get_QO()
                            ->get_mul_nucomp()
                            ->get_RelativeGenerator()
                            ->template to_log<S>();

  update_distance_add(Distance, Distance, relative_distance_1);
  update_distance_add(Distance, Distance, relative_distance_2);
}

// START: TEMPORARY SECTION FOR REGULATORLENSTRADATA METHODS
template <class T, class S>
void QuadraticInfElement<T, S>::adjust(const ZZ &a) {

    S bound = to<S>(a);

    if (Distance > bound) {
      while (Distance > bound) {
        inverse_rho();
      }
    }

    else {
      while (Distance <= bound) {
        baby_step();
      }
      inverse_rho();
    }
    return;
}

template <class T, class S>
void QuadraticInfElement<T, S>::assign(const HashEntryReal<T>) {}

template <class T, class S> void QuadraticInfElement<T, S>::assign_one() {}

template <class T, class S>
QuadraticInfElement<T, S> QuadraticInfElement<T, S>::conjugate() const {}

template <class T, class S> ZZ QuadraticInfElement<T, S>::eval() {}

template <class T, class S>
HashEntryReal<T> QuadraticInfElement<T, S>::hash_real() const {}

template <class T, class S> bool QuadraticInfElement<T, S>::is_one() {}

template <class T, class S> void QuadraticInfElement<T, S>::inverse_rho() {
  T a = qib.get_a(), b = qib.get_b(), c = qib.get_c(), q, a2, r, temp, nb, oa;
  S relative_generator, relative_distance;

    oa = a;
    a = -c;
    a2 = a << 1;

    // q = floor((rootD + b) / 2a)
    temp = FloorRootDelta + b;
    DivRem(q, r, temp, a2);
    if (temp < 0 && !IsZero(r)) {
      --q;
    }

  // Computation of the relative distance
  relative_generator = abs(to<S>(2 * oa) / (to<S>(b) - sqrt(to<S>(Delta))));
  relative_distance = log(relative_generator);

  // Updating the distance
  update_distance_subtract(Distance, Distance, relative_distance);

    // b = 2aq - b
    nb = FloorRootDelta - r;

    // c =  q*((nb - b)/2) - c
    c = q * ((nb - b) >> 1) - oa;

    b = nb;
    if (a < 0) {
      a = -a;
      c = -c;
    }

    qib.assign(a, b, c);
}

template <class T, class S>
S QuadraticInfElement<T, S>::get_baby_steps(
    IndexedHashTable<HashEntryReal<T>> &prin_list, const ZZ &B,
    const QuadraticInfElement<T, S> &A) {}

template <class T, class S>
S QuadraticInfElement<T, S>::get_baby_steps(
    IndexedHashTable<HashEntryReal<T>> &prin_list, const ZZ &B,
    const QuadraticInfElement<T, S> &A, long l, long &M){

};

//   qo_distance<T>
//   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
//                  const ZZ &B, const qi_pair<T> &A) {};
//   qo_distance<T>
//   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
//                  const ZZ &B, const qi_pair<T> &A, long l, long &M) {};j

} // namespace ANTL
