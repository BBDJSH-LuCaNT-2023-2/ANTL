#ifndef REGULATOR_LENSTRA_DATA_ZZ_H
#define REGULATOR_LENSTRA_DATA_ZZ_H

#include <ANTL/Quadratic/Regulator/RegulatorLenstra.hpp>

using namespace NTL;
using namespace ANTL;

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class U> class RegulatorLenstraData<ZZ, U> {
private:

  //DBG_CONSTANTS
  bool DBG_LENSTR = false;
  bool DBG_EHRERR = false;
  bool DBG_GOQCNM = false;
  bool DBG_BSGSGL = false;
  bool DBG_BSGSES = false;
  bool DBG_IPLIST = false;
  bool DBG_SHANKS = false;
  bool DBG_APPRHR = false;
  bool DBG_OPTIMK = false;

  U regulator;

  QuadraticOrder<ZZ> *quadratic_order;

  ZZ delta;

  L_function<ZZ> *l_function;

  bool parallel;

  const long OQvals_cnum[20] = {2269,   5741,   10427,  16183,  22901,
                                30631,  39209,  48731,  59063,  70237,
                                82223,  95009,  108571, 122921, 137983,
                                153817, 170341, 187631, 205589, 224261};

  const long OQvals[20] = {947,   2269,  3929,  6011,  8447,  11093, 14149,
                           17393, 20921, 24733, 28807, 33151, 37619, 42533,
                           47507, 52859, 58321, 64231, 70099, 76463};

  ZZ max_memory{4000000000};

  IndexedHashTable<HashEntryReal<ZZ, U>> prin_list;
  long prinlist_l; // distance between consecutive table entries

  long prinlist_M; // max # of baby-steps between table entries

  ZZ prinlist_s; // giant-step distance

  bool Rbsgs; // true if R was computed using BSGS

  bool Rconditional; // true if correctness of R relies on ERH

public:
  RegulatorLenstraData(QuadraticOrder<ZZ> *quadratic_order_arg,
                       L_function<ZZ> *l_function_arg);

  void regulator_lenstra();

  U get_regulator();

  ZZ lower_bound_hR();

private:
  ZZ estimate_hR_error();

  long get_optimal_Q_cnum();

  long bsgs_getl(const ZZ &K, ZZ &N, ZZ &entry_size, RR &mu, bool nodist);

  RR get_mu(const ZZ &delta);

  void bsgs_getentrysize(ZZ &entry_size, bool nodist);

  void init_prinlist(const ZZ &N, long l, ZZ &s, long &M,
                     QuadraticInfElement<ZZ, U> &G);

  void regulator_bsgs(ZZ &bound);

  ZZ approximate_hR();

  void optimize_K(ZZ &bound, const ZZ &S, const ZZ &N, long l);

  void combine_BSGS(U &dist, const QuadraticInfElement<ZZ, U> &DD,
                    const HashEntryReal<ZZ, U> *F);

  void combine_conj_BSGS(U &dist, const QuadraticInfElement<ZZ, U> &DD,
                         const HashEntryReal<ZZ, U> *F);

  long get_optimal_Q();

  long generate_optimal_Q();

};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class U>
RegulatorLenstraData<ZZ, U>::RegulatorLenstraData(
    QuadraticOrder<ZZ> *quadratic_order_arg, L_function<ZZ> *l_function_arg) {

  quadratic_order = quadratic_order_arg;
  delta = quadratic_order->get_discriminant();

  l_function = l_function_arg;

  parallel = false;
  Rbsgs = false;        // true if R was computed using BSGS
  Rconditional = false; // true if correctness of R relies on ERH
}

// RegulatorLenstraData<ZZ, U>::regulator_lenstra
// Task: Computes the regulator using an O(D^1/5) baby-step giant-step algorithm
// of Shanks and improvements of Lenstra.

template <class U>
void RegulatorLenstraData<ZZ, U>::regulator_lenstra() {

  //
  // initialize hash table
  //
  ZZ K, N, B, entry_size, u, s, s2;
  long l = 1, M = 1;
  RR mu;
  QuadraticInfElement<ZZ, U> AA{*quadratic_order}, G{*quadratic_order};

  U S;
  // qo_distance<ZZ> S;
  clear(S);

  if(DBG_LENSTR) std::cout << "LENSTR: into EHRERR" << std::endl;
  K = estimate_hR_error() >> 1;
  if(DBG_LENSTR) std::cout << "LENSTR: out of EHRERR" << std::endl;
  if(DBG_LENSTR) std::cout << "LENSTR: K is " << K << std::endl;
  if(DBG_LENSTR) std::cout << "LENSTR: into BSGSGL" << std::endl;
  l = bsgs_getl(K, N, entry_size, mu, false);
  if(DBG_LENSTR) std::cout << "LENSTR: out of BSGSGL" << std::endl;

  if(DBG_LENSTR) std::cout << "LENSTR: into IPLIST" << std::endl;
  init_prinlist(N, l, s, M, G);
  if(DBG_LENSTR) std::cout << "LENSTR: out of IPLIST" << std::endl;

  B = N * l;
  if (IsZero(B)) {
    if(DBG_LENSTR) std::cout << "LENSTR: into SHANKS" << std::endl;
    regulator_bsgs(B);
    if(DBG_LENSTR) std::cout << "LENSTR: out of SHANKS" << std::endl;
    S = regulator;
  }

  if(DBG_LENSTR) std::cout << "LENSTR: regulator is " << regulator << std::endl;
  if (IsZero(S)) {
    //
    // compute approximation of hR
    //
    if(DBG_LENSTR) std::cout << "LENSTR: into APPRHR" << std::endl;
    ZZ E = approximate_hR();
    if(DBG_LENSTR) std::cout << "LENSTR: out of APPRHR" << std::endl;
    nuclose(AA, E);

    if (AA.is_one())
      S = AA.get_distance();
  }
  //
  // compute list of baby steps (distance < B)
  //
  QuadraticInfElement<ZZ, U> A, C, CC, D, DD, GG;
  HashEntryReal<ZZ, U> *F;

  if (IsZero(S)) {

    A.assign_one();
    if (l == 1)
      regulator = G.get_baby_steps(prin_list, B, A);
    else
      regulator = G.get_baby_steps(prin_list, B, A, l, M);

    if (!IsZero(regulator)) {
      nuclose(C, FloorToZZ(regulator));
      regulator = C.get_distance();
      Rbsgs = true;
    }

    prinlist_M = M;

    if (IsZero(regulator)) {
      G.adjust(s);

      u = (s << 1);
      G.giant_step(G);
      //sqr(G, G); makeshift square above
      G.adjust(u);
      if (G.is_one())
        G.baby_step();

      conjugate(GG, G);
      while (abs(GG.eval()) <= u)
        GG.inverse_rho();
      while (abs(GG.eval()) > u)
        GG.baby_step();

      s2 = s = AA.eval();
      C = AA;
      D = AA;
    }
  }
  if(DBG_LENSTR) std::cout << "LENSTR: regulator is " << regulator << std::endl;
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

      C.giant_step(G);
      //mul(C, C, G); makeshift multiplication above
      C.adjust(s);

      s2 -= u;

      D.giant_step(GG);
      //mul(D, D, GG); makeshift multiplication above
      D.adjust(s2);
    }
  }
  if(DBG_LENSTR) std::cout << "LENSTR: regulator is " << regulator << std::endl;

  //
  // factor out h*
  //

  if (IsZero(regulator)) {

    // verify that Regulator > S^2/3 (= B) using BS-GS
    ZZ B, N, entry_size;
    RR mu;
    long l;
    RR temp = log(to_RR(S)) * to_RR (2) / to_RR (3);

    conv(B, ceil(exp(temp)));

    l = bsgs_getl(B, N, entry_size, mu, true);
    optimize_K(B, FloorToZZ(S), N, l);

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

    nuclose(C, FloorToZZ(regulator));
    regulator = C.get_distance();
  }
}

template <class U> U RegulatorLenstraData<ZZ, U>::get_regulator() {
  return regulator;
}

// RegulatorLenstraData<ZZ, U>::estimate_hR_error(RegulatorLenstraData<ZZ, U>
// &rl_data)
// Task: returns L such that |hR - hR'| < exp(L)^2

template <class U> ZZ RegulatorLenstraData<ZZ, U>::estimate_hR_error() {
  ZZ err;
  RR Aval, Fval, temp;

  if (l_function->terms_used(1) == 0){
    if(DBG_EHRERR) std::cout << "EHRERR: terms_used == 0 !" << std::endl;
    return ZZ::zero();
  }

  long n = get_optimal_Q_cnum();
  RR FI = l_function->approximateL1(n);

  if(DBG_EHRERR) std::cout << "EHRERR: n is " << n << std::endl;
  if(DBG_EHRERR) std::cout << "EHRERR: FI is " << FI << std::endl;

  Aval = l_function->calculate_L1_error(delta, l_function->terms_used(1));
  Fval = exp(Aval) - 1;

  if(DBG_EHRERR) std::cout << "EHRERR: Aval is " << Aval << std::endl;
  if(DBG_EHRERR) std::cout << "EHRERR: Fval is " << Fval << std::endl;

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

  if(DBG_EHRERR) std::cout << "EHRERR: Fval is " << Fval << std::endl;

  err = CeilToZZ(Fval * log(to_RR(2))) >> 1;
  return err;
}

template <class U> long RegulatorLenstraData<ZZ, U>::get_optimal_Q_cnum() {
  long dlog;
  ZZ temp;

  temp = FloorToZZ(NTL::log10(to_RR(abs(delta))));
  conv(dlog, temp);

  return OQvals_cnum[(dlog / 5)];
}

template <class U>
long RegulatorLenstraData<ZZ, U>::bsgs_getl(const ZZ &K, ZZ &N, ZZ &entry_size,
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
  if (DBG_BSGSGL) std::cout << "BSGSGL: K is "  << K << std::endl;
  if (DBG_BSGSGL) std::cout << "BSGSGL: mu is " << mu << std::endl;
  if (DBG_BSGSGL) std::cout << "BSGSGL: n is "  << n << std::endl;
  temp = floor(SqrRoot(to_RR(K) * mu / (to_RR(1) * n)));
  conv(N, temp);
  l = 1;
  if (DBG_BSGSGL) std::cout << "BSGSGL: N is " << N << std::endl;
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
  if (DBG_BSGSGL) std::cout << "BSGSGL: N is " << N << std::endl;
  return l;
}

// quadratic_order<ZZ>::get_mu()
// Task: determines the number of baby-steps which can be done in the time of
// one giant step.  These values are read in from a table (file) which is
// computed at compile-time.

template <class U> RR RegulatorLenstraData<ZZ, U>::get_mu(const ZZ &delta) {
  //  return to_RR(2)*((to_RR(NumBits(delta) - 5) / to_RR(10)) + to_RR(6));
  //  return (to_RR(NumBits(delta) - 5) / to_RR(10)) + to_RR(6);
  return to_RR(6.85) + to_RR(10.62) * to_RR(NumBits(delta)) / to_RR(135);
}

template <class U>
void RegulatorLenstraData<ZZ, U>::bsgs_getentrysize(ZZ &entry_size,
                                                    bool nodist) {
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

// quadratic_order<ZZ>::init_prinlist
// Task: initializes the list of reduced principal ideals for BS-GS
// computations. Uses current contents if appropriate.

template <class U>
void RegulatorLenstraData<ZZ, U>::init_prinlist(const ZZ &N, long l, ZZ &s,
                                               long &M,
                                               QuadraticInfElement<ZZ, U> &G) {
  long lsize, P;
  if(DBG_IPLIST) std::cout << "IPLIST: Begin" << std::endl;

  // compute B and s
  if(DBG_IPLIST) std::cout << "IPLIST: N is " << N << std::endl;
  if(DBG_IPLIST) std::cout << "IPLIST: l is " << l << std::endl;

  P = 3 + NumBits(SqrRoot(delta));
  s = N * l - P;

  if(DBG_IPLIST) std::cout << "IPLIST: P is " << P << std::endl;
  if(DBG_IPLIST) std::cout << "IPLIST: s is " << s << std::endl;
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
    if(DBG_IPLIST) std::cout << "IPLIST! 1.1" << std::endl;
    conv(lsize, N + P);
    if(DBG_IPLIST) std::cout << "IPLIST! 1.2" << std::endl;
    lsize += 100;
    if(DBG_IPLIST) std::cout << "IPLIST! 1.3" << std::endl;
    prin_list.initialize(lsize);
    if(DBG_IPLIST) std::cout << "IPLIST! 1.4" << std::endl;
    G.assign_one();
    if(DBG_IPLIST) std::cout << "IPLIST! 1.5" << std::endl;
    prinlist_l = l;
    if(DBG_IPLIST) std::cout << "IPLIST! 1.6" << std::endl;
    prinlist_s = s;
    if(DBG_IPLIST) std::cout << "IPLIST! 1.7" << std::endl;
  }
  if(DBG_IPLIST) std::cout << "IPLIST: s is " << s << std::endl;
  if(DBG_IPLIST) std::cout << "IPLIST: End" << std::endl;
}

//
// quadratic_order<ZZ>::regulator_bsgs
// Task: Computes the regulator using an O(sqrt(R)) baby-step giant-step
// algorithm.
//

template <class U>
void RegulatorLenstraData<ZZ, U>::regulator_bsgs(ZZ &bound) {

  // initialize hash table
  ZZ K, B, N, entry_size, u, s;
  RR mu;
  long l, M = 1;

  QuadraticInfElement<ZZ, U> G{*quadratic_order};

//   if(DBG_SHANKS) std::cout << "SHANKS: Begin" << std::endl;
//   if(DBG_SHANKS) std::cout << "K is " << K << std::endl;
//   if(DBG_SHANKS) std::cout << "delta is " << delta << std::endl;

  if (bound == 0)
    // replacing get_bound for now...
    K = SqrRoot(delta);
  // get_bound(K);
  else
    K = bound;

//   if(DBG_SHANKS) std::cout << "K is " << K << std::endl;

//   if(DBG_SHANKS) std::cout << "SHANKS: into BSGSGL" << std::endl;
  l = bsgs_getl(K, N, entry_size, mu, false);
//   if(DBG_SHANKS) std::cout << "SHANKS: out of BSGSGL" << std::endl;

//   if(DBG_SHANKS) std::cout << "SHANKS: into IPLIST" << std::endl;
  init_prinlist(N, l, s, M, G);
//   if(DBG_SHANKS) std::cout << "SHANKS: out of IPLIST" << std::endl;

  B = N * l;
//   if(DBG_SHANKS) std::cout << "SHANKS: Begin baby step list computation" << std::endl;
  // compute list of baby steps (distance < B)
  QuadraticInfElement<ZZ, U> A{*quadratic_order}, C{*quadratic_order}, AA{*quadratic_order};
  HashEntryReal<ZZ, U> *F;

  if(DBG_SHANKS) std::cout << "SHANKS: regulator is " << regulator << std::endl;

  A.assign_one();
//   if(DBG_SHANKS) std::cout << "SHANKS: Begin get_baby_steps" << std::endl;
//   if(DBG_SHANKS) std::cout << "SHANKS: l == " << l << std::endl;
  if (l == 1)
    regulator = G.get_baby_steps(prin_list, B, A);
  else
    regulator = G.get_baby_steps(prin_list, B, A, l, M);

//   if(DBG_SHANKS) std::cout << "SHANKS: End get_baby_steps" << std::endl;
  prinlist_M = M;
  if(DBG_SHANKS) std::cout << "SHANKS: regulator is " << regulator << std::endl;

  if (!IsZero(regulator)) {
    Rbsgs = true;
  }

  if (IsZero(regulator)) {
    Rbsgs = false;
//     if(DBG_SHANKS) std::cout << "SHANKS: Begin adjust" << std::endl;
    G.adjust(s);
//     if(DBG_SHANKS) std::cout << "SHANKS: End adjust" << std::endl;
    u = (s << 1);

    A = G;
//     if(DBG_SHANKS) std::cout << "SHANKS: Begin giant_step" << std::endl;
    if(DBG_SHANKS) std::cout << "SHANKS: (" << G.get_qib().get_a() << ", "
                    << G.get_qib().get_b() << ", "
                    << G.get_qib().get_c() << ") "
                    << G.get_distance() << std::endl;

    G.giant_step(G);

    if(DBG_SHANKS) std::cout << "SHANKS: (" << G.get_qib().get_a() << ", "
                    << G.get_qib().get_b() << ", "
                    << G.get_qib().get_c() << ") "
                    << G.get_distance() << std::endl;
//     if(DBG_SHANKS) std::cout << "SHANKS: End giant_step" << std::endl;
    //sqr(G, G); makeshift square above
//     if(DBG_SHANKS) std::cout << "SHANKS: Begin adjust, baby_step" << std::endl;
    if(DBG_SHANKS) std::cout << "SHANKS: G distance is " << G.get_distance() << std::endl;
    G.adjust(u);
    if(DBG_SHANKS) std::cout << "SHANKS: G distance is " << G.get_distance() << std::endl;
    if (G.is_one())
      G.baby_step();
//     if(DBG_SHANKS) std::cout << "SHANKS: End adjust, baby_step" << std::endl;
    if(DBG_SHANKS) std::cout << "SHANKS: G distance is " << G.get_distance() << std::endl;
  }

//   if(DBG_SHANKS) std::cout << "SHANKS: End baby step list computation" << std::endl;
  //
  // compute giant steps until R is found or the bound is exceeded
  //
//   if(DBG_SHANKS) std::cout << "SHANKS: Begin giant step traversal" << std::endl;
  long i;

//   if(DBG_SHANKS) std::cout << "SHANKS: s and u are " << s << " and " << u << std::endl;
  while (IsZero(regulator) && (bound == 0 || A.eval() <= bound)) {
    s += u;


    if(DBG_SHANKS) std::cout << "SHANKS: (" << G.get_qib().get_a() << ", "
                    << G.get_qib().get_b() << ", "
                    << G.get_qib().get_c() << ") "
                    << G.get_distance() << std::endl;

    if(DBG_SHANKS) std::cout << "SHANKS: (" << A.get_qib().get_a() << ", "
                    << A.get_qib().get_b() << ", "
                    << A.get_qib().get_c() << ") "
                    << A.get_distance() << std::endl;
    A.giant_step(G);

    if(DBG_SHANKS) std::cout << "SHANKS: (" << A.get_qib().get_a() << ", "
                    << A.get_qib().get_b() << ", "
                    << A.get_qib().get_c() << ") "
                    << A.get_distance() << std::endl;

                    //mul(A, A, G); makeshift multiplication above
    A.adjust(s);

    // search for A, rho_1(A), ..., rho_l(A) in the hash table
    if(DBG_SHANKS) std::cout << "SHANKS: A distance is " << A.get_distance() << std::endl;
    AA = A;
    if(DBG_SHANKS) std::cout << "SHANKS: AA distance is " << AA.get_distance() << std::endl;
    if(DBG_SHANKS) std::cout << "SHANKS: M  is " << M << std::endl;
    if(DBG_SHANKS) std::cout << "SHANKS: regulator is " << regulator << std::endl;
    for (i = 0; i < M && IsZero(regulator); ++i) {
      F = prin_list.search(AA.hash_real());
      if (F) {
        // found AA in the hash table!
        if(DBG_SHANKS) std::cout << "SHANKS: found AA in the hash table! " << std::endl;
        combine_BSGS(regulator, AA, F);
        if(DBG_SHANKS) std::cout << "SHANKS: regulator is " << regulator << std::endl;
        nuclose(C, FloorToZZ(regulator));
        C.adjust(regulator);
        regulator = C.get_distance();
        if(DBG_SHANKS) std::cout << "SHANKS: regulator is " << regulator << std::endl;

      } else {
        F = prin_list.search((AA.conjugate()).hash_real());
        if (F) {
          // found A^-1 in the hash table!
          if(DBG_SHANKS) std::cout << "SHANKS: found AA^-1 in the hash table! " << std::endl;
          combine_conj_BSGS(regulator, AA, F);
          nuclose(C, FloorToZZ(regulator));
          C.adjust(regulator);
          regulator = C.get_distance();
        } else if (i < M)
          // AA.rho(); Swapping this out for AA.baby_step();
          if(DBG_SHANKS) std::cout << "SHANKS: AA distance is " << AA.get_distance() << std::endl;
          AA.baby_step();
          if(DBG_SHANKS) std::cout << "SHANKS: AA distance is " << AA.get_distance() << std::endl;
      }
    }
    if(DBG_SHANKS) std::cout << "SHANKS: regulator is " << regulator << std::endl;

  }
  if(DBG_SHANKS) std::cout << "SHANKS: End giant step traversal" << std::endl;
  if (!IsZero(regulator))
    Rconditional = false;

  if(DBG_SHANKS) std::cout << "SHANKS: End" << std::endl;
}

//
// quadratic_order<ZZ>::approximate_hR()
//
// Task:
//      returns an approximation of hR
//

template <class U> ZZ RegulatorLenstraData<ZZ, U>::approximate_hR() {
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

template <class U>
void RegulatorLenstraData<ZZ, U>::optimize_K(ZZ &bound, const ZZ &S, const ZZ &N,
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

template <class U>
void RegulatorLenstraData<ZZ, U>::combine_BSGS(
    U &dist, const QuadraticInfElement<ZZ, U> &DD, const HashEntryReal<ZZ, U> *F) {
  dist = DD.get_distance() - F->get_d();
}

template <class U>
void RegulatorLenstraData<ZZ, U>::combine_conj_BSGS(
    U &dist, const QuadraticInfElement<ZZ, U> &DD, const HashEntryReal<ZZ, U> *F) {
  dist = DD.get_distance() + F->get_d() - deg(DD.get_qib().get_a());
}


//
// lower_bound_hR()
// Task: returns a lower bound of hR such that L < hR < 2L
//
template <class U>
ZZ RegulatorLenstraData<ZZ, U>::lower_bound_hR() {
  ZZ hR;
  RR temp, FI;

  //START: Temporary variables needed (previously declared in ANTL-Import's quadratic_order
  bool unconditional = true;
  int info = 0;
  bool use_tables = false;
  //FINISH: Temporary variables needed (previously declared in ANTL-Import's quadratic_order

  if (unconditional) {
    if (info > 1) {
      cout << "Lower bound hR with L(0,X)" << endl;
    }

    // compute hR apporximation

    if (quadratic_order->is_imaginary()) {
      if (use_tables)
        FI = l_function->approximateL0_ImaginaryNumberField_table(
            log(sqrt(double(2))));
      else
        FI = l_function->approximateL0_ImaginaryNumberField(log(sqrt(double(2))));

      if (info > 1) {
        cout << "L(0,X) approx " << FI << endl;
      }

      temp = FI;
      if (delta == -4)
        temp *= 2;
      if (delta == -3)
        temp *= 3;

    } else {
      if (use_tables)
        FI = l_function->approximateL0_RealNumberField_table(log(sqrt(double(2))));
      else
        FI = l_function->approximateL0_RealNumberField(log(sqrt(double(2))));

      if (info > 1) {
        cout << "L(0,X) approx " << FI << endl;
      }
      temp = FI;
    }
  } else {
    if (info > 1) {
      cout << "Lower bound hR with L(1,X)" << endl;
    }

    if (use_tables)
      FI = l_function->approximateL1_table();
    else {
      long n = get_optimal_Q();
      if (info > 2)
        cout << "using n = " << n << endl;
      FI = l_function->approximateL1(n);
    }

    if (info > 1) {
      cout << "L(1,X) approx " << FI << endl;
    }

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
  }

  hR = CeilToZZ(temp / SqrRoot(to_RR(2)));

  if (info > 1)
    cout << "Lower bound = " << hR << endl;

  return hR;
}

//
// quadratic_order<ZZ>::get_optimal_Q()
//
// Task: returns a value of Q from a pre-computed table which will compute h*
// such that h* < h < 2h*.
//

template <class U>
long RegulatorLenstraData<ZZ, U>::get_optimal_Q() {
  long Dlog;
  ZZ temp;

  temp = FloorToZZ(log10(to_RR(abs(delta))));
  conv(Dlog, temp);

  if ((Dlog / 5) > 19)
    return generate_optimal_Q();

  return OQvals[Dlog / 5];
}

//
// quadratic_order<ZZ>::generate_optimal_Q()
// Task: computes the optimal value of Q for computing h* such that h* < h < 2h*
//

template <class U> long RegulatorLenstraData<ZZ, U>::generate_optimal_Q() {
  long OQ;
  RR A, l2;
  PrimeSeq primes;

  l2 = log(sqrt(to_RR(2)));

  // set OQ = 3
  OQ = primes.next();

  do {
    OQ = primes.next();
    A = l_function->calculate_L1_error(delta, OQ);
  } while (A >= l2);

  return OQ;
}
#endif
