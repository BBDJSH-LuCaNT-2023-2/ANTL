#ifndef REGULATOR_LENSTRA_DATA_H
#define REGULATOR_LENSTRA_DATA_H

#include <ANTL/common.hpp>

#include <ANTL/HashTable/HashEntryReal.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>
#include <ANTL/L_function/L_function.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>
#include <ANTL/Quadratic/QuadraticOrder.hpp>

NTL_CLIENT
using namespace ANTL;

namespace ANTL {
template <class T> class RegulatorLenstraData {
private:
  RR regulator;

  QuadraticOrder<T> *quadratic_order;

  T delta;

  L_function<T> *l_function;

  bool parallel;

  const long OQvals_cnum[20] = {2269,   5741,   10427,  16183,  22901,
                                30631,  39209,  48731,  59063,  70237,
                                82223,  95009,  108571, 122921, 137983,
                                153817, 170341, 187631, 205589, 224261};

  ZZ max_memory{4000000000};

  IndexedHashTable<HashEntryReal<T>> prin_list;

  long prinlist_l; // distance between consecutive table entries

  long prinlist_M; // max # of baby-steps between table entries

  ZZ prinlist_s; // giant-step distance

  bool Rbsgs; // true if R was computed using BSGS

  bool Rconditional; // true if correctness of R relies on ERH

public:
  RegulatorLenstraData(QuadraticOrder<ZZ> *quadratic_order_arg,
                       L_function<ZZ> *l_function_arg);

  void regulator_lenstra();

private:
  ZZ estimate_hR_error();

  long get_optimal_Q_cnum();

  long bsgs_getl(const ZZ &K, ZZ &N, ZZ &entry_size, RR &mu, bool nodist);

  RR get_mu(const T &delta);

  void bsgs_getentrysize(ZZ &entry_size, bool nodist);

  void init_prinlist(const ZZ &N, long l, ZZ &s, long &M,
                     QuadraticInfElement<T> &G);

  void regulator_bsgs(ZZ &bound);

  ZZ approximate_hR();

  void optimize_K(ZZ &bound, const ZZ &S, const ZZ &N, long l);

  RR func(const RR &K, const RR &S, const RR &N, const RR &G, const RR &n,
          long l);

  RR dfunc(const RR &K, const RR &S, const RR &N, const RR &G, const RR &n,
           long l);

  void nuclose(QuadraticInfElement<ZZ> &C, const ZZ &n);

  void combine_BSGS(RR &dist, const QuadraticInfElement<T> &DD,
                    const HashEntryReal<T> *F);

  void combine_conj_BSGS(RR &dist, const QuadraticInfElement<T> &DD,
                         const HashEntryReal<T> *F);
};

// START: Forward Declaring template specializations

template <>
RegulatorLenstraData<ZZ>::RegulatorLenstraData(
    QuadraticOrder<ZZ> *quadratic_order_arg, L_function<ZZ> *l_function_arg);

template <> ZZ RegulatorLenstraData<ZZ>::estimate_hR_error();

template <> long RegulatorLenstraData<ZZ>::get_optimal_Q_cnum();

template <> RR RegulatorLenstraData<ZZ>::get_mu(const ZZ &delta);

template <>
void RegulatorLenstraData<ZZ>::bsgs_getentrysize(ZZ &entry_size, bool nodist);

template <> ZZ RegulatorLenstraData<ZZ>::approximate_hR();

template <>
void RegulatorLenstraData<ZZ>::nuclose(QuadraticInfElement<ZZ> &C, const ZZ &n);

// FINISH:Forward Declaring template specializations

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <>
RegulatorLenstraData<ZZ>::RegulatorLenstraData(
    QuadraticOrder<ZZ> *quadratic_order_arg, L_function<ZZ> *l_function_arg) {

  quadratic_order = quadratic_order_arg;
  delta = quadratic_order->get_discriminant();

  l_function = l_function_arg;

  parallel = false;
  Rbsgs = false;        // true if R was computed using BSGS
  Rconditional = false; // true if correctness of R relies on ERH
}

// RegulatorLenstraData<T>::regulator_lenstra
// Task: Computes the regulator using an O(D^1/5) baby-step giant-step algorithm
// of Shanks and improvements of Lenstra.

template <class T> void RegulatorLenstraData<T>::regulator_lenstra() {

  //
  // initialize hash table
  //

  ZZ K, N, B, entry_size, u, s, s2;
  long l = 1, M = 1;
  RR mu;
  QuadraticInfElement<T> AA, G;

  RR S;
  // qo_distance<T> S;
  clear(S);

  K = estimate_hR_error() >> 1;
  l = bsgs_getl(K, N, entry_size, mu, false);
  init_prinlist(N, l, s, M, G);
  B = N * l;
  if (IsZero(B)) {
    regulator_bsgs(B);
    S = regulator;
  }

  if (IsZero(S)) {
    //
    // compute approximation of hR
    //

    ZZ E = approximate_hR();
    nuclose(AA, E);

    if (AA.is_one())
      S = AA.get_distance();
  }

  //
  // compute list of baby steps (distance < B)
  //

  QuadraticInfElement<T> A, C, CC, D, DD, GG;
  HashEntryReal<T> *F;

  if (IsZero(S)) {

    A.assign_one();
    if (l == 1)
      regulator = G.get_baby_steps(prin_list, B, A);
    else
      regulator = G.get_baby_steps(prin_list, B, A, l, M);

    if (!IsZero(regulator)) {
      nuclose(C, FloorToZZ(log(regulator) / log(RR(2))));
      regulator = C.get_distance();
      Rbsgs = true;
    }

    prinlist_M = M;

    if (IsZero(regulator)) {
      G.adjust(s);

      u = (s << 1);
      sqr(G, G);
      G.adjust(u);
      if (G.is_one())
        G.baby_step();

      GG = G.conjugate();
      while (abs(GG.eval()) <= u)
        GG.inverse_rho();
      while (abs(GG.eval()) > u)
        GG.baby_step();

      s2 = s = AA.eval();
      C = AA;
      D = AA;
    }
  }

  //
  // compute giant steps until R is found or the bound is exceeded
  //

  long i;

  while (IsZero(regulator) && IsZero(S)) {

    // search for C, rho_1(C), ..., rho_l(C) in the hash table
    CC = C;
    DD = D;
    for (i = 0; i < M && IsZero(S); ++i) {
      F = prin_list.search(CC.hash_real());
      if (F) {
        // found CC in the hash table!

        combine_BSGS(S, CC, F);
      } else {
        F = prin_list.search((CC.conjugate()).hash_real());
        if (F) {
          // found CC^-1 in the hash table!

          combine_conj_BSGS(S, CC, F);
        } else {
          F = prin_list.search(DD.hash_real());
          if (F) {
            // found DD in the hash table!

            combine_BSGS(S, DD, F);
          } else {
            F = prin_list.search((DD.conjugate()).hash_real());
            if (F) {
              // found DD^-1 in the hash table!

              combine_conj_BSGS(S, DD, F);
            } else if (i < M) {
              CC.baby_step();
              DD.baby_step();
            }
          }
        }
      }
    }

    if (IsZero(S)) {
      s += u;
      mul(C, C, G);
      C.adjust(s);

      s2 -= u;
      mul(D, D, GG);
      D.adjust(s2);
    }
  }

  //
  // factor out h*
  //

  if (IsZero(regulator)) {

    // verify that Regulator > S^2/3 (= B) using BS-GS
    ZZ B, N, entry_size;
    RR mu;
    long l;
    RR temp = log(to_RR(FloorToZZ(log(S) / log(RR(2))))) * to_RR(2) / to_RR(3);

    conv(B, ceil(exp(temp)));

    l = bsgs_getl(B, N, entry_size, mu, true);
    optimize_K(B, FloorToZZ(log(to_RR(FloorToZZ(log(S) / log(RR(2)))))), N, l);

    regulator_bsgs(B);

    ZZ hstar, Pmax;

    //     if (info > 1 && !IsZero(regulator)) {
    //       hstar = (FloorToZZ(log(S) / log(RR(2)))) /
    //       (FloorToZZ(log(regulator) / log(RR(2))));
    //
    //       cout << "Found R = " << regulator << endl;
    //       cout << "h* = " << hstar << endl;
    //     } else {
    //       Pmax = 1 + FloorToZZ(log(S) / log(RR(2))) / RR(B);
    //       find_hstar(hstar, S, Pmax, t1);
    //     }

    nuclose(C, FloorToZZ(log(regulator) / log(RR(2))));
    regulator = C.get_distance();
  }
}

// RegulatorLenstraData<ZZ>::estimate_hR_error(RegulatorLenstraData<ZZ>
// &rl_data)
// Task: returns L such that |hR - hR'| < exp(L)^2

template <> ZZ RegulatorLenstraData<ZZ>::estimate_hR_error() {
  ZZ err;
  RR Aval, Fval, temp;

  if (l_function->terms_used(1) == 0)
    return ZZ::zero();

  long n = get_optimal_Q_cnum();
  RR FI = l_function->approximateL1(n);

  Aval = l_function->calculate_L1_error(delta, l_function->terms_used(1));
  Fval = exp(Aval) - 1;
  temp = 1 - exp(-Aval);
  if (temp > Fval)
    Fval = temp;

  if (quadratic_order->is_imaginary()) {
    // h = sqrt(delta) L / Pi
    Fval *= FI * SqrRoot(to_RR(-delta)) / ComputePi_RR();
    if (delta == -4)
      temp *= 2;
    if (delta == -3)
      temp *= 3;
  } else {
    // hR = sqrt(delta) L / 2
    Fval *= FI * SqrRoot(to_RR(delta)) / 2;
  }

  err = CeilToZZ(Fval * log(to_RR(2))) >> 1;
  return err;
}

template <> long RegulatorLenstraData<ZZ>::get_optimal_Q_cnum() {
  long dlog;
  ZZ temp;

  temp = FloorToZZ(log10(to_RR(abs(quadratic_order->get_discriminant()))));
  conv(dlog, temp);

  return OQvals_cnum[(dlog / 5)];
}

template <class T>
long RegulatorLenstraData<T>::bsgs_getl(const ZZ &K, ZZ &N, ZZ &entry_size,
                                        RR &mu, bool nodist) {
  ZZ maxN, rootK;
  RR n, temp;
  long l;

  // get mu (# baby steps per giant step)
  mu = get_mu(delta);

  // n = number of slave processes
  if (parallel)
    n = parallel;
  else
    set(n);

  // compute max number of baby steps to store
  bsgs_getentrysize(entry_size, nodist);
  maxN = 1 + (max_memory / (to_ZZ(n) * entry_size));

  // compute number of baby steps assuming l=1, N=sqrt(KG/2n)
  temp = floor(SqrRoot(to_RR(K) * mu / (to_RR(1) * n)));
  conv(N, temp);
  l = 1;

  if (N > maxN) {
    // compute l = sqrt(KG/2n) / N;
    //    N = maxN;
    entry_size *= 3;
    entry_size >>= 1;
    N = 1 + (max_memory / (to_ZZ(n) * entry_size));
    temp = floor((SqrRoot(to_RR(K) * mu / (to_RR(2) * n)) / to_RR(N)));
    conv(l, temp);
    if (l < 1)
      l = 1;
  }

  return l;
}

// quadratic_order<ZZ>::get_mu()
// Task: determines the number of baby-steps which can be done in the time of
// one giant step.  These values are read in from a table (file) which is
// computed at compile-time.

template <> RR RegulatorLenstraData<ZZ>::get_mu(const ZZ &delta) {
  //  return to_RR(2)*((to_RR(NumBits(delta) - 5) / to_RR(10)) + to_RR(6));
  //  return (to_RR(NumBits(delta) - 5) / to_RR(10)) + to_RR(6);
  return to_RR(6.85) + to_RR(10.62) * to_RR(NumBits(delta)) / to_RR(135);
}

template <>
void RegulatorLenstraData<ZZ>::bsgs_getentrysize(ZZ &entry_size, bool nodist) {
  ZZ temp;

  entry_size = 3 * NTL_BITS_PER_LONG; // space for pointers
  if (nodist) {
    entry_size += SqrRoot(delta).size() * NTL_ZZ_NBITS; // a coeff
    entry_size += 8;                                    // b coeff
  } else {
    // temp = to_ZZ(1) << qo_distance<ZZ>::get_p();
    temp = to_ZZ(1) << 64;
    entry_size +=
        temp.size() * NTL_ZZ_NBITS + 3 * SqrRoot(delta).size() * NTL_ZZ_NBITS;
  }

  // convert bits to bytes
  entry_size >>= 2;
}

// quadratic_order<T>::init_prinlist
// Task: initializes the list of reduced principal ideals for BS-GS
// computations. Uses current contents if appropriate.

template <class T>
void RegulatorLenstraData<T>::init_prinlist(const ZZ &N, long l, ZZ &s, long &M,
                                            QuadraticInfElement<T> &G) {
  long lsize, P;

  // compute B and s
  P = 3 + NumBits(SqrRoot(delta));
  s = N * l - P;

  // initialize hash table
  if (prin_list.no_of_elements() > 0) {
    G.assign(prin_list.last_entry());
    l = prinlist_l;
    M = prinlist_M;
    if (prinlist_s > s)
      s = prinlist_s;
    else
      prinlist_s = s;
  } else {
    conv(lsize, N + P);
    lsize += 100;
    prin_list.initialize(lsize);
    G.assign_one();
    prinlist_l = l;
    prinlist_s = s;
  }
}

//
// quadratic_order<T>::regulator_bsgs
// Task: Computes the regulator using an O(sqrt(R)) baby-step giant-step
// algorithm.
//

template <class T> void RegulatorLenstraData<T>::regulator_bsgs(ZZ &bound) {

  // initialize hash table
  ZZ K, B, N, entry_size, u, s;
  RR mu;
  long l, M = 1;

  QuadraticInfElement<T> G;

  if (bound == 0)
    // replacing get_bound for now...
    bound = SqrRoot(delta);
  // get_bound(K);
  else
    K = bound;

  l = bsgs_getl(K, N, entry_size, mu, false);
  init_prinlist(N, l, s, M, G);
  B = N * l;

  // compute list of baby steps (distance < B)
  QuadraticInfElement<T> A, C, AA;
  HashEntryReal<T> *F;

  A.assign_one();
  if (l == 1)
    regulator = G.get_baby_steps(prin_list, B, A);
  else
    regulator = G.get_baby_steps(prin_list, B, A, l, M);

  prinlist_M = M;

  if (!IsZero(regulator)) {
    Rbsgs = true;
  }

  if (IsZero(regulator)) {
    Rbsgs = false;

    G.adjust(s);

    u = (s << 1);
    A = G;
    sqr(G, G);
    G.adjust(u);
    if (G.is_one())
      G.baby_step();
  }

  //
  // compute giant steps until R is found or the bound is exceeded
  //

  long i;

  while (IsZero(regulator) && (bound == 0 || A.eval() <= bound)) {
    s += u;
    mul(A, A, G);
    A.adjust(s);

    // search for A, rho_1(A), ..., rho_l(A) in the hash table
    AA = A;
    for (i = 0; i < M && IsZero(regulator); ++i) {
      F = prin_list.search(AA.hash_real());
      if (F) {
        // found AA in the hash table!

        combine_BSGS(regulator, AA, F);

        nuclose(C, FloorToZZ(log(regulator) / log(RR(2))));
        regulator = C.get_distance();

      } else {
        F = prin_list.search((AA.conjugate()).hash_real());
        if (F) {
          // found A^-1 in the hash table!

          combine_conj_BSGS(regulator, AA, F);
          nuclose(C, FloorToZZ(log(regulator) / log(RR(2))));
          regulator = C.get_distance();
        } else if (i < M)
          // AA.rho(); Swapping this out for AA.baby_step();
          AA.baby_step();
      }
    }
  }

  if (!IsZero(regulator))
    Rconditional = false;
}

//
// quadratic_order<ZZ>::approximate_hR()
//
// Task:
//      returns an approximation of hR
//

template <> ZZ RegulatorLenstraData<ZZ>::approximate_hR() {
  ZZ hR;
  RR temp, FI;

  long n = get_optimal_Q_cnum();

  FI = l_function->approximateL1(n);

  if (quadratic_order->is_imaginary()) {
    // h = sqrt(delta) L / Pi
    temp = FI * SqrRoot(to_RR(-delta)) / ComputePi_RR();
    if (delta == -4)
      temp *= 2;
    if (delta == -3)
      temp *= 3;
  } else {
    // hR = sqrt(delta) L / 2
    temp = FI * SqrRoot(to_RR(delta)) / 2;
  }

  hR = CeilToZZ(temp);

  return hR;
}

template <class T>
void RegulatorLenstraData<T>::optimize_K(ZZ &bound, const ZZ &S, const ZZ &N,
                                         long l) {
  RR K = to_RR(bound);
  RR mu = get_mu(delta);
  RR rS = to_RR(S);
  RR rN = to_RR(N);
  RR rn;

  if (parallel)
    rn = to_RR(parallel);
  else
    set(rn);

  RR F, dF;

  F = func(K, rS, rN, mu, rn, l);
  dF = dfunc(K, rS, rN, mu, rn, l);
  cout << "K = " << K << ", F(K) = " << F << ", F'(K) = " << dF << endl;

  while (abs(F) > to_RR(0.00001)) {
    K = K - F / dF;
    F = func(K, rS, rN, mu, rn, l);
    dF = dfunc(K, rS, rN, mu, rn, l);
    //    cout << "K = " << K << ", F(K) = " << F << ", F'(K) = " << dF << endl;
  }
  cout << "K = " << K << ", F(K) = " << F << endl;

  bound = FloorToZZ(K);
}

inline RR func(const RR &K, const RR &S, const RR &N, const RR &G, const RR &n,
               long l) {
  RR val;

  //  cout << "\nF:  K = " << K << ", S = " << S << endl;

  //  val = (1.5*(G+1)*S*log(K)) / (n*K*(log(S)-log(K)));
  val = ((G + 1) * S * (log(K) - log(N * l))) / (n * K * (log(S) - log(K)));
  val -= SqrRoot((2 * K * G) / n);
  if (l > 1)
    val -= K / (2 * n * N);

  return val;
}

inline RR dfunc(const RR &K, const RR &S, const RR &N, const RR &G, const RR &n,
                long l) {
  RR val, temp;

  temp = log(S) - log(K);

  //  val = (1.5*(G+1)*S*(1-log(K))) / (K*K*n*temp);
  //  val += (1.5*(G+1)*S*log(K)) / (K*K*n*temp*temp);
  val = ((G + 1) * S * (1 - log(K) + log(N * l))) / (K * K * n * temp);
  val += ((G + 1) * S * (log(K) - log(N * l))) / (K * K * n * temp * temp);
  val -= sqrt(G / (2 * K * n));
  if (l > 1)
    val -= to_RR(1) / (2 * N * n);
  return val;
}

template <>
void RegulatorLenstraData<ZZ>::nuclose(QuadraticInfElement<ZZ> &C,
                                       const ZZ &n) {
  long i, k = 0;
  ZZ j, ex, s;

  C.assign_one();
  if (IsZero(n))
    return;

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
    sqr(C, C);

    if (IsOdd(j))
      ++s;

    C.adjust(s);

    j >>= 1;
  }
}

template <class T>
void RegulatorLenstraData<T>::combine_BSGS(RR &dist,
                                           const QuadraticInfElement<T> &DD,
                                           const HashEntryReal<T> *F) {
  dist = DD.get_distance() - F->get_d();
}

template <class T>
void RegulatorLenstraData<T>::combine_conj_BSGS(
    RR &dist, const QuadraticInfElement<T> &DD, const HashEntryReal<T> *F) {
  dist = DD.get_distance() + F->get_d() - deg(DD.get_qib().get_a());
}

} // namespace ANTL
#endif
