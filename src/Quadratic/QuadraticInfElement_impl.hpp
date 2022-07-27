using namespace ANTL;

namespace ANTL {

template <class T, class S> QuadraticInfElement<T, S>::QuadraticInfElement(){};

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

  // Computation of the QuadraticIdealBase via a single step in the continued
  // fraction expansion
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
                            ->get_mul_comp()
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
void QuadraticInfElement<T, S>::adjust(const S &bound) {

  S current_difference, previous_difference;
  current_difference = bound - Distance;
//   std::cout << "ADJUST taraget distance is " << bound << std::endl;
//   std::cout << qib << " with distance" << Distance << std::endl;
  if (-0.000001 >= bound - Distance) {
    while (Distance > bound) {
      inverse_rho();
      previous_difference = current_difference;
      current_difference = bound - Distance;
//       std::cout << qib << " with distance" << Distance << std::endl;
    }
    if (abs(previous_difference) < abs(current_difference)) {
      baby_step();
    }
      return;
  }

  else {
    // Testing Floating point equivalence is hard; a simple Distance <= bound
    // won't work reliably. This workaround ensures correctness within a
    // reasonable (read arbritaty) accuracy level
    while (-0.000001 < bound - Distance) {
      baby_step();
      previous_difference = current_difference;
      current_difference = bound - Distance;
//       std::cout << qib << " with distance" << Distance << std::endl;
    }
    inverse_rho();
    previous_difference = current_difference;
    current_difference = bound - Distance;
    if (abs(previous_difference) < abs(current_difference)) {
      baby_step();
    }
//     std::cout << qib << " with distance" << Distance << std::endl;
  }


  return;
}
/*
template <class T, class S>
void QuadraticInfElement<T, S>::adjust_to_one() {

  QuadraticInfElement<T, S> forward_qie{*qib.get_QO()}, backward_qie{*qib.get_QO()};

  forward_qie.qib = qib;
  forward_qie.Distance = Distance;

  backward_qie.qib = qib;
  backward_qie.Distance = Distance;

  while(!(forward_qie.qib.IsOne() || backward_qie.qib.IsOne())) {
    forward_qie.baby_step();
    backward_qie.inverse_rho();
  }

  if(forward_qie.qib.IsOne()) {
    qib = forward_qie.get_qib();
    Distance = forward_qie.Distance;
  } else {
    qib = backward_qie.get_qib();
    Distance = backward_qie.Distance;
  }

  return;
}*/

//
// qi_pair<T>::assign(quadratic_form<T>)
// Task: set to a copy of B.
//

template <class T, class S>
void QuadraticInfElement<T, S>::assign(const QuadraticIdealBase<T> &B) {
  qib.assign(B);
  qib.reduce();
  clear(Distance);
}

template <class T, class S>
void QuadraticInfElement<T, S>::assign(const HashEntryReal<T, S> &her_a) {
  T a = her_a.get_a();
  T b = her_a.get_b();
  T c = (b * b - Delta) / (a << 2);

  qib.assign(a, b, c);
  qib.normalize();

  // this->reduce_real ();
  Distance = her_a.get_d();
}

template <class T, class S> void QuadraticInfElement<T, S>::assign_one() {
  qib.assign_one();
  Distance = 0;
}

template <class T, class S> ZZ QuadraticInfElement<T, S>::eval() {
  return FloorToZZ(Distance);
}

template <class T, class S>
HashEntryReal<T, S> QuadraticInfElement<T, S>::hash_real() const {
  // below is the implementatio for class T = ZZ
  // it is here in the generic implementation temporarily
  // eventually a QuadraticInfElement<ZZ, S> class definition will be required

  return HashEntryReal<T, S>(get_qib().get_a(), get_qib().get_b(), Distance);
}

template <class T, class S> bool QuadraticInfElement<T, S>::is_one() const {
  return qib.IsOne();
}

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

//
// get_baby_steps (qo_hash_entry_real)
//
// Adds all baby steps (starting with *this) with distance <= B to the hash
// table
//
template <class T, class S>
S QuadraticInfElement<T, S>::get_baby_steps(
    IndexedHashTable<HashEntryReal<T, S>> &prin_list, const ZZ &B,
    const QuadraticInfElement<T, S> &A) {
  // below is the implementation for class T = ZZ
  // it is here in the generic implementation temporarily
  // eventually a QuadraticInfElement<ZZ, S> class definition will be required
  ZZ old_a = qib.get_a(), old_b = qib.get_b(), old_c = qib.get_c();
  ZZ q, r, new_c, new_a, new_b;
  S relative_generator, relative_distance, regulator, initial_distance = Distance;

  prin_list.hash(this->hash_real());

  do {
    // Computation of the QuadraticIdealBase via a single step in the continued
    // fraction expansion
    DivRem(q, r, old_b + FloorRootDelta, 2 * old_a);

    new_c = -old_a;
    new_b = FloorRootDelta - r;
    new_a = q * ((old_b - new_b) / 2) - old_c;

    // Computation of the relative distance
    relative_generator =
        abs(to<S>(2*new_a) / (to<S>(new_b) - sqrt(to<S>(Delta))));
    relative_distance = log(relative_generator);

    // The two checks below are for the cases where consecutive coefficients are
    // found to be equal in the continued fraction expansion
    if (old_a == new_a) {
      regulator = 2 * (Distance);
//       std::cout << "old_a == new_a" << std::endl;
      update_distance_add(regulator, regulator, relative_distance);
      regulator = regulator - log(to<S>(old_a));
      return regulator;
    }

    if (old_b == new_b && (!IsOne(old_a))) {
//       std::cout << "old_b == new_b" << std::endl;
      regulator = 2 * (Distance);
      regulator = regulator - log(to<S>(old_a));
      return regulator;
    }

    qib.assign(new_a, new_b, new_c);
    qib.normalize();

    // Updating the distance
    update_distance_add(Distance, Distance, relative_distance);

    old_a = new_a;
    old_b = new_b;
    old_c = new_c;

    if (qib == A.get_qib()) {
      return Distance;
    }

    prin_list.hash(hash_real());

  } while (-0.000000001 <= to<double>(B) - Distance);

  regulator = 0;
  return regulator;
}

//
// get_baby_steps (l, qo_hash_entry_real)
//
// Adds all baby steps (starting with *this) with distance <= B to the hash
// table. Only those ideals close to multiples of l are added to the hash table.
// M is the maximum number of baby-steps between two ideals added to the hash
// table.
//
template <class T, class S>
S QuadraticInfElement<T, S>::get_baby_steps(
    IndexedHashTable<HashEntryReal<T, S>> &prin_list, const ZZ &B,
    const QuadraticInfElement<T, S> &A, long l, long &M) {
  // below is the implementatio for class T = ZZ
  // it is here in the generic implementation temporarily
  // eventually a QuadraticInfElement<ZZ, S> class definition will be required

  ZZ old_a = qib.get_a(), old_b = qib.get_b(), old_c = qib.get_c();
  ZZ q, r, new_c, new_a, new_b;

  ZZ temp, a2, nb, na;

  S relative_generator, relative_distance, regulator, new_d;
  ZZ sl;
  long s, currM;

  ZZ bound = to_ZZ(10000000);

  sl = FloorToZZ(Distance) + l;
  currM = 0;

  new_d = Distance;

  do {

    // Computation of the QuadraticIdealBase via a single step in the continued
    // fraction expansion
    DivRem(q, r, old_b + FloorRootDelta, 2 * old_a);

    new_c = -old_a;
    new_b = FloorRootDelta - r;
    new_a = q * ((old_b - new_b) / 2) - old_c;

    update_distance_multiply(new_d, new_d, to<S>(new_b));
    new_d = new_d / to<S>(2 * old_a);

    if (Distance <= sl && new_d > sl) {
      if (currM > M)
        M = currM;
      currM = 0;

      if (sl > bound) {
        cout << "bound = " << bound << endl;
        bound += 10000000;
      }

      do {
        prin_list.hash(this->hash_real());
        sl += l;
      } while (sl < new_d);
    } else
      ++currM;

    // The two checks below are for the cases where consecutive coefficients are
    // found to be equal in the continued fraction expansion
    if (A.is_one() && old_a == new_a) {
      regulator = 2 * Distance;
      update_distance_add(regulator, Distance, relative_distance);
      regulator = regulator / to<S>(old_a);
      return regulator;
    }

    if (A.is_one() && old_b == new_b && (!IsOne(old_a))) {
//       std::cout << "old_b == new_b" << std::endl;
      regulator = 2 * Distance;
      regulator = regulator / to<S>(old_a);
      return regulator;
    }

    qib.assign(new_a, new_b, new_c);
    qib.normalize();

    // Updating the distance
    update_distance_add(Distance, Distance, relative_distance);

    old_a = new_a;
    old_b = new_b;
    old_c = new_c;

    if (qib == A.get_qib())
      return Distance;

  } while (Distance <= B);

  regulator = 0;
  return regulator;
}

//   qo_distance<T>
//   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
//                  const ZZ &B, const qi_pair<T> &A) {};
//   qo_distance<T>
//   get_baby_steps(indexed_hash_table<qo_hash_entry_real<T>> &prin_list,
//                  const ZZ &B, const qi_pair<T> &A, long l, long &M) {};

template <class T, class S>
void conjugate(QuadraticInfElement<T, S> &qie_a,
               QuadraticInfElement<T, S> const &qie_b) {

  qie_a = qie_b;
  qie_a.qib.set_b(-qie_b.qib.get_b());

  update_distance_add(qie_a.Distance, qie_a.get_distance(), qie_b.get_distance());
  update_distance_negate(qie_a.Distance, qie_a.Distance);
  update_distance_subtract(qie_a.Distance, qie_a.Distance,
                           log(to<S>(qie_a.qib.get_a())));

  qie_a.qib.normalize();
}

template <class T, class S>
QuadraticInfElement<T, S> QuadraticInfElement<T, S>::conjugate() {

  QuadraticInfElement<T, S> qie_conj{*qib.get_QO()};

  qie_conj.qib = qib;
  qie_conj.qib.set_b(-qib.get_b());


  update_distance_add(qie_conj.Distance, Distance, qie_conj.Distance);
  update_distance_negate(qie_conj.Distance, qie_conj.Distance);
  update_distance_subtract(qie_conj.Distance, qie_conj.Distance,
                           log(to<S>(qie_conj.qib.get_a())));

  qie_conj.qib.normalize();

  return qie_conj;
}

template <class T, class S>
bool operator==(const QuadraticInfElement<T, S> &qie_a,
                const QuadraticInfElement<T, S> &qie_b) {
  return qie_a.qib == qie_a.qib;
}

template <class T, class S>
void nuclose(QuadraticInfElement<T, S> &C, const ZZ &n) {
  long i, k = 0;
  ZZ j, ex, s;
  C.assign_one();

  if (IsZero(n)) {
    return;
  }
  // compute binary expansion of ex (hi order to low order)
  ex = abs(n);
  clear(j);
  while (!IsOne(ex)) {
    j <<= 1;
    if (IsOdd(ex))
      ++j;
    ex >>= 1;
    ++k;
  }

  s = 1;
  C.adjust(s);


  for (i = 1; i <= k; ++i) {
    s <<= 1;
    C.giant_step(C);

    // sqr(C, C); makeshift square above

    if (IsOdd(j))
      ++s;
    C.adjust(s);

    j >>= 1;
  }
}

} // namespace ANTL
