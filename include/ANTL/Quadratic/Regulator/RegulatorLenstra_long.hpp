#ifndef REGULATOR_LENSTRA_DATA_LONG_H
#define REGULATOR_LENSTRA_DATA_LONG_H

#include <list>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra.hpp>

using namespace NTL;
using namespace ANTL;

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class U> class RegulatorLenstraData<long, U> {
private:
  // DBG_CONSTANTS
  bool DBG_LENSTR = false;
  bool DBG_EHRERR = false;
  bool DBG_GOQCNM = false;
  bool DBG_BSGSGL = false;
  bool DBG_GETMU_ = false;
  bool DBG_BSGSES = false;
  bool DBG_IPLIST = false;
  bool DBG_SHANKS = false;
  bool DBG_APPRHR = false;
  bool DBG_OPTIMK = false;
  bool DBG_LOBOHR = false;
  bool DBG_GETOPQ = false;
  bool DBG_GENOPQ = false;
  bool DBG_FHSTAR = false;

// public:
//   void CHG_DBG_LENSTR(bool flag) {DBG_LENSTR = flag;}
//   void CHG_DBG_EHRERR(bool flag) {DBG_EHRERR = flag;}
//   void CHG_DBG_GOQCNM(bool flag) {DBG_GOQCNM = flag;}
//   void CHG_DBG_BSGSGL(bool flag) {DBG_BSGSGL = flag;}
//   void CHG_DBG_GETMU_(bool flag) {DBG_GETMU_ = flag;}
//   void CHG_DBG_BSGSES(bool flag) {DBG_BSGSES = flag;}
//   void CHG_DBG_IPLIST(bool flag) {DBG_IPLIST = flag;}
//   void CHG_DBG_SHANKS(bool flag) {DBG_SHANKS = flag;}
//   void CHG_DBG_APPRHR(bool flag) {DBG_APPRHR = flag;}
//   void CHG_DBG_OPTIMK(bool flag) {DBG_OPTIMK = flag;}
//   void CHG_DBG_LOBOHR(bool flag) {DBG_LOBOHR = flag;}
//   void CHG_DBG_GETOPQ(bool flag) {DBG_GETOPQ = flag;}
//   void CHG_DBG_GENOPQ(bool flag) {DBG_GENOPQ = flag;}
//   void CHG_DBG_FHSTAR(bool flag) {DBG_FHSTAR = flag;}
//
// private:
  U regulator;

  ZZ hstar;

  QuadraticOrder<long> *quadratic_order;

  long delta;

  L_function<long> *l_function;

  bool parallel;

  const long OQvals_cnum[20] = {2269,   5741,   10427,  16183,  22901,
                                30631,  39209,  48731,  59063,  70237,
                                82223,  95009,  108571, 122921, 137983,
                                153817, 170341, 187631, 205589, 224261};

  const long OQvals[20] = {947,   2269,  3929,  6011,  8447,  11093, 14149,
                           17393, 20921, 24733, 28807, 33151, 37619, 42533,
                           47507, 52859, 58321, 64231, 70099, 76463};

  ZZ max_memory{8000000000};

  IndexedHashTable<HashEntryReal<long, U>> prin_list;

  long prinlist_l; // distance between consecutive table entries

  long prinlist_M; // max # of baby-steps between table entries

  U prinlist_s; // giant-step distance

  bool Rbsgs; // true if R was computed using BSGS

  bool Rconditional; // true if correctness of R relies on ERH

  std::string case_type = "";

public:
  RegulatorLenstraData(QuadraticOrder<long> *quadratic_order_arg,
                       L_function<long> *l_function_arg);

  void regulator_lenstra();

  void regulator_bsgs(ZZ &bound);

  U get_regulator();

  ZZ get_hstar() {return hstar;}

  RR lower_bound_hR();

  void set_case_type(std::string found_case_type);

  std::string get_case_type();

  U approximate_hR();

private:
  ZZ estimate_hR_error();

  long get_optimal_Q_cnum();

  long bsgs_getl(const ZZ &K, ZZ &N, ZZ &entry_size, RR &mu, bool nodist);

  RR get_mu(const long &delta);

  void bsgs_getentrysize(ZZ &entry_size, bool nodist);

  void init_prinlist(const ZZ &N, long l, U &s, long &M,
                     QuadraticInfElement<long, U> &G);

  void optimize_K(ZZ &bound, const ZZ &S, const ZZ &N, long l);

  void combine_BSGS(U &dist, const QuadraticInfElement<long, U> &DD,
                    const HashEntryReal<long, U> *F);

  void combine_conj_BSGS(U &dist, const QuadraticInfElement<long, U> &DD,
                         const HashEntryReal<long, U> *F);

  long get_optimal_Q();

  long generate_optimal_Q();

  void find_hstar(ZZ &hstar, const U &S, const ZZ &Pmax);
};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class U>
RegulatorLenstraData<long, U>::RegulatorLenstraData(
    QuadraticOrder<long> *quadratic_order_arg, L_function<long> *l_function_arg) {

  quadratic_order = quadratic_order_arg;
  delta = quadratic_order->get_discriminant();

  l_function = l_function_arg;

  parallel = false;
  Rbsgs = false;        // true if R was computed using BSGS
  Rconditional = false; // true if correctness of R relies on ERH
}

template <class U>
void RegulatorLenstraData<long, U>::set_case_type(std::string found_case_type) {
  case_type += found_case_type;
}

template <class U> std::string RegulatorLenstraData<long, U>::get_case_type() {
  return case_type;
}

// RegulatorLenstraData<long, U>::regulator_lenstra
// Task: Computes the regulator using an O(D^1/5) baby-step giant-step algorithm
// of Shanks and improvements of Lenstra.

template <class U> void RegulatorLenstraData<long, U>::regulator_lenstra() {
  std::cout << "reg_len<long> is using reg_len<long>.hpp" << std::endl;

  //
  // initialize hash table
  //
  ZZ K, N, B, entry_size;
  U u, s, s2;
  long l = 1, M = 1;
  RR mu;
  QuadraticInfElement<long, U> AA{*quadratic_order}, G{*quadratic_order};

  U S;
  // qo_distance<ZZ> S;
  clear(S);
  clear(regulator);

  U E = approximate_hR();

  K = estimate_hR_error() >> 1;

  if (DBG_LENSTR) {
    std::cout << "LENSTR: K is " << K << std::endl;
  }

  // l is baby-step increment for which are stored in the table
  l = bsgs_getl(K, N, entry_size, mu, false);

  init_prinlist(N, l, s, M, G);

  if (DBG_LENSTR) {
    std::cout << "LENSTR: s is " << s << std::endl;
    std::cout << "LENSTR: N is " << N << std::endl;
    std::cout << "LENSTR: l is " << l << std::endl;
  }

  B = N * l;
  if (DBG_LENSTR) {
    std::cout << "LENSTR: B is " << B << std::endl;
  }

  if (IsZero(B)) {
    regulator_bsgs(B);
    S = regulator;
    if (!IsZero(regulator)) {
      set_case_type(" -- Initial regulator_bsfs -- ");
    }
  }

  if (DBG_LENSTR) {
    std::cout << "LENSTR: regulator is " << regulator << std::endl;
    std::cout << "LENSTR: S is " << S << std::endl;
    std::cout << "LENSTR: E is " << E << std::endl;
  }

  if (IsZero(S)) {
    //
    // compute approximation of hR
    //

    nuclose(AA, FloorToZZ(E));
    AA.adjust(E);

    if (AA.is_one()) {
      S = AA.get_distance();
      if (DBG_LENSTR) {
        std::cout << "LENSTR: AA.is_one() - setting S!" << std::endl;
        std::cout << "LENSTR: S is " << S << std::endl;
      }
    }
  }

  //
  // compute list of baby steps (distance < B)
  //
  QuadraticInfElement<long, U> A{*quadratic_order}, C{*quadratic_order},
      C_gap_check{*quadratic_order}, D{*quadratic_order},
      D_gap_check{*quadratic_order}, GG{*quadratic_order}, C1{*quadratic_order},
      D1{*quadratic_order};
  HashEntryReal<long, U> *F;

  if (IsZero(S)) {

    A.assign_one();
    if (l == 1) {
      regulator = G.get_baby_steps(prin_list, B, A);
      if (DBG_LENSTR) {
        std::cout << "LENSTR: G.get_baby_steps(prin_list, B, A) " << std::endl;
        std::cout << "LENSTR: B is " << B << std::endl;
      }
    } else {
      regulator = G.get_baby_steps(prin_list, B, A, l, M);
      if (DBG_LENSTR) {
        std::cout << "LENSTR: G.get_baby_steps(prin_list, B, A, l, M) "
                  << std::endl;
      }
    }

    if (DBG_LENSTR) {
        std::cout << "LENSTR: prin_list is " << prin_list << std::endl;
    }

    if (!IsZero(regulator)) {
      Rbsgs = true;
    }

    prinlist_M = M;

    if (IsZero(regulator)) {
      G.adjust(s);
      if (DBG_LENSTR) {
        std::cout << "LENSTR: G.adjust(s) " << std::endl;
        std::cout << "LENSTR: s is " << s << std::endl;
      }
      u = 2 * s;
      if (DBG_LENSTR) {
        std::cout << "LENSTR: u is " << u << std::endl;
      }
      G.giant_step(G);
      if (DBG_LENSTR) {
        std::cout << "LENSTR: G.giant_step(G) " << std::endl;
      }
      // sqr(G, G); makeshift square above
      G.adjust(u);
      if (DBG_LENSTR) {
        std::cout << "LENSTR: G.adjust(s) " << std::endl;
      }
      if (G.is_one()) {
        G.baby_step();
        if (DBG_LENSTR) {
          std::cout << "LENSTR: G.baby_step() " << std::endl;
        }
      }

      conjugate(GG, G);
      if (DBG_LENSTR) {
        std::cout << "LENSTR: conjugate(GG, G) " << std::endl;
      }
      // return;
      while (abs(GG.eval()) <= u) {
        GG.inverse_rho();
        if (DBG_LENSTR) {
          std::cout << "LENSTR: GG.inverse_rho() " << std::endl;
        }
      }
      while (abs(GG.eval()) > u) {
        if (DBG_LENSTR) {
          std::cout << "LENSTR: GG.eval()) > u" << std::endl;
        }

        if (DBG_LENSTR) {
          std::cout << "LENSTR: GG is (" << GG.get_qib().get_a() << ", "
                    << GG.get_qib().get_b() << " , " << GG.get_qib().get_c()
                    << ") with distance " << GG.get_distance() << std::endl;
          std::cout << "LENSTR: u is " << u << std::endl;
        }

        GG.baby_step();
        if (DBG_LENSTR) {
          std::cout << "LENSTR: GG.baby_step() " << std::endl;
        }
      }

      if (DBG_LENSTR) {
        std::cout << "AA is " << AA.get_qib() << " with distance "
                  << AA.get_distance() << std::endl;
      }
      s2 = s = AA.get_distance();
      C = AA;
      D = AA.conjugate();
    }

    else {
      set_case_type(" -- baby_step list construction -- ");
    }
  }

  if (DBG_LENSTR) {
    std::cout << "LENSTR: regulator is " << regulator << std::endl;
    std::cout << "LENSTR: S is " << S << std::endl;
  }

  //
  // compute giant steps until R is found or the bound is exceeded
  //

  // Saving initial ideals
  C1 = C;
  D1 = D;

  // Perform giant-steps on both C and D until one is found in the
  // baby-step list
  while (IsZero(regulator) && IsZero(S)) {

    if (DBG_LENSTR) {
      std::cout << "C distance is " << C.get_distance() << std::endl;
      std::cout << "D distance is " << D.get_distance() << std::endl;
    }

    C_gap_check = C;
    D_gap_check = D;

    // If there are gaps of M baby-steps in our baby-step list, we must perform
    // M baby-steps on C and D and check for each of these in the list
    for (long i = 0; i < M && IsZero(S); ++i) {

      F = prin_list.search(C_gap_check.hash_real());
      if (F) {

        // found C in the hash table!
        set_case_type(" -- C found -- ");
        if (DBG_LENSTR) {
          std::cout << "found C in the hash table!" << std::endl;
        }

        combine_BSGS(S, C_gap_check, F);
        break;
      }

      F = prin_list.search(D_gap_check.hash_real());
      if (F) {

        // found D in the hash table!
        set_case_type(" -- D found -- ");
        if (DBG_LENSTR) {
          std::cout << "found D in the hash table!" << std::endl;
        }

        // combine_conj_BSGS(S, D_gap_check, F);
        S = AA.get_distance() -
            ((D_gap_check.get_distance() - D1.get_distance()) - F->get_d()) -
            log(ZZ(AA.get_qib().get_a()));
        break;
      }

      C_gap_check.baby_step();
      D_gap_check.baby_step();
    }

    if(abs(S) < 1) {
      S = 0;
    }

    // If nothing has been found, S will still be zero; Continue with giant
    // steps
    if (IsZero(S)) {
      s += u;
      s2 -= u;

      // mul(C, C, G); makeshift multiplication below
      C.giant_step(G);
      C.adjust(s);

      // mul(D, D, GG); makeshift multiplication below
      D.giant_step(G);
      // D.adjust(s2);
    }
  }

  if (DBG_LENSTR) {
    std::cout << "LENSTR: regulator is " << regulator << std::endl;
    std::cout << "LENSTR: S is " << S << std::endl;
  }

  //
  // factor out h*
  //
  if (S < 0) {
    S *= -1;
  }
  if (IsZero(regulator)) {

    if (DBG_LENSTR) {
      std::cout << "LENSTR: Factoring out h* " << std::endl;
      std::cout << "LENSTR: regulator is " << regulator << std::endl;
      std::cout << "LENSTR: S is " << S << std::endl;
    }

    // verify that Regulator > S^2/3 (= B) using BS-GS
    ZZ B, N, entry_size;
    RR mu;
    long l;
    if (DBG_LENSTR) {
      std::cout << "LENSTR: Before log(S), S is " << S << std::endl;
    }
    RR temp = log(to_RR(S)) * to_RR(2) / to_RR(3);

    if (DBG_LENSTR) {
      std::cout << "LENSTR: After log(S)" << std::endl;
      std::cout << "LENSTR: Before conv(B, ceil(exp(temp)))" << std::endl;
    }

    conv(B, ceil(exp(temp)));
    if (DBG_LENSTR) {
      std::cout << "LENSTR: After conv(B, ceil(exp(temp)))" << std::endl;
    }

    l = bsgs_getl(B, N, entry_size, mu, true);

    if (DBG_LENSTR) {
      std::cout << "LENSTR: B is " << B << std::endl;
      std::cout << "LENSTR: FloorToZZ(S) is " << FloorToZZ(S) << std::endl;
      std::cout << "LENSTR: N is " << N << std::endl;
      std::cout << "LENSTR: l is " << l << std::endl;
    }

    // THE BELOW LINE SHOULD EVENTUALLY BE UNCOMMENTED
    //     optimize_K(B, FloorToZZ(S), N, l);

    if (DBG_LENSTR) {
      std::cout << "LENSTR: After optimize_K(B, FloorToZZ(S), N, l)"
                << std::endl;
      std::cout << "LENSTR: B is " << B << std::endl;
    }

    // THE BELOW LINE SHOULD EVENTUALLY BE UNCOMMENTED
    B = CeilToZZ(E / SqrRoot(to_RR(K)));
    regulator_bsgs(B);
    if(!IsZero(regulator)){
      set_case_type(" -- second regulator_bsgs -- ");
    }

    if (DBG_LENSTR) {
      std::cout << "LENSTR: After regulator_bsgs(B) regulator is " << regulator
                << std::endl;
      std::cout << "LENSTR: S is " << S << std::endl;
    }

    ZZ Pmax;
    if (!IsZero(regulator)) {
      hstar = (FloorToZZ(log(S) / log(RR(2)))) /
              (FloorToZZ(log(regulator) / log(RR(2))));
      if (DBG_LENSTR) {
        std::cout << "Found R = " << regulator << std::endl;
        std::cout << "h* = " << hstar << std::endl;
      }
    } else {
      Pmax = 1 + CeilToZZ(to_RR(S) / to_RR(B));
//             Pmax = FloorToZZ(sqrt(S));
      if (DBG_LENSTR) {
        std::cout << "LENSTR: Pmax is " << Pmax << std::endl;
      }
      find_hstar(hstar, S, Pmax);
    }

    if (DBG_LENSTR) {
      std::cout << "LENSTR: Factoring out h* - S is " << S << std::endl;
      std::cout << "LENSTR: Factoring out h* - hstar is " << hstar << std::endl;
      std::cout << "LENSTR: Factoring out h* - regulator is " << regulator
                << std::endl;
      std::cout << "LENSTR: Factoring out h* - nuclose start" << std::endl;
    }

    nuclose(C, FloorToZZ(regulator));
    C.adjust(regulator);

    if (DBG_LENSTR) {
      std::cout << "LENSTR: Factoring out h* - nuclose finish" << std::endl;
      std::cout << "LENSTR: Factoring out h* - assign regulator" << std::endl;
    }

    regulator = C.get_distance();
  }

  if (DBG_LENSTR) {
    std::cout << "LENSTR: Factoring out h* - END S is" << S << std::endl;
    std::cout << "LENSTR: Factoring out h* - END regulator is" << regulator
              << std::endl;
  }
}

template <class U> U RegulatorLenstraData<long, U>::get_regulator() {
  return regulator;
}

// RegulatorLenstraData<long, U>::estimate_hR_error(RegulatorLenstraData<long, U>
// &rl_data)
// Task: returns L such that |hR - hR'| < exp(L)^2

template <class U> ZZ RegulatorLenstraData<long, U>::estimate_hR_error() {
  if (DBG_LENSTR || DBG_EHRERR) {
    std::cout << "EHRERR: START" << std::endl;
  }
  ZZ err;
  RR Aval, Fval, temp;

  if (DBG_EHRERR)
    std::cout << "EHRERR: l_function->terms_used(1) is "
              << l_function->terms_used(1) << std::endl;

  long n = get_optimal_Q_cnum();
  RR FI = l_function->approximateL1(n);

  if (l_function->terms_used(1) == 0) {
    if (DBG_EHRERR)
      std::cout << "EHRERR: terms_used == 0 !" << std::endl;
    return ZZ::zero();
  }

  if (DBG_EHRERR)
    std::cout << "EHRERR: n is " << n << std::endl;
  if (DBG_EHRERR)
    std::cout << "EHRERR: FI is " << FI << std::endl;

  Aval = l_function->calculate_L1_error(delta, l_function->terms_used(1));
  Fval = exp(Aval) - 1;

  if (DBG_EHRERR)
    std::cout << "EHRERR: Aval is " << Aval << std::endl;
  if (DBG_EHRERR)
    std::cout << "EHRERR: Fval is " << Fval << std::endl;

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

  if (DBG_EHRERR)
    std::cout << "EHRERR: Fval is " << Fval << std::endl;

  err = CeilToZZ(Fval * log(to_RR(2))) >> 1;

  if (DBG_LENSTR || DBG_EHRERR) {
    std::cout << "EHRERR: FINISH" << std::endl;
  }
  return err;
}

template <class U> long RegulatorLenstraData<long, U>::get_optimal_Q_cnum() {
  if (DBG_LENSTR || DBG_GOQCNM) {
    std::cout << "DBG_GOQCNM: START" << std::endl;
  }

  long dlog;
  ZZ temp;

  temp = FloorToZZ(NTL::log10(to_RR(abs(delta))));
  conv(dlog, temp);

  if (DBG_LENSTR || DBG_GOQCNM) {
    std::cout << "DBG_GOQCNM: FINISH" << std::endl;
  }
  return OQvals_cnum[(dlog / 5)];
}

template <class U>
long RegulatorLenstraData<long, U>::bsgs_getl(const ZZ &K, ZZ &N, ZZ &entry_size,
                                            RR &mu, bool nodist) {
  if (DBG_LENSTR || DBG_BSGSGL) {
    std::cout << "BSGSGL: START" << std::endl;
  }
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
  if (DBG_BSGSGL)
    std::cout << "BSGSGL: K is " << K << std::endl;
  if (DBG_BSGSGL)
    std::cout << "BSGSGL: mu is " << mu << std::endl;
  if (DBG_BSGSGL)
    std::cout << "BSGSGL: n is " << n << std::endl;
  temp = floor(SqrRoot(to_RR(K) * mu / (to_RR(1) * n)));
  conv(N, temp);
  l = 1;
  if (DBG_BSGSGL)
    std::cout << "BSGSGL: N is " << N << std::endl;
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
  if (DBG_BSGSGL)
    std::cout << "BSGSGL: N is " << N << std::endl;
  if (DBG_LENSTR || DBG_BSGSGL) {
    std::cout << "BSGSGL: FINISH" << std::endl;
  }
  return l;
}

// quadratic_order<ZZ>::get_mu()
// Task: determines the number of baby-steps which can be done in the time of
// one giant step.  These values are read in from a table (file) which is
// computed at compile-time.

template <class U> RR RegulatorLenstraData<long, U>::get_mu(const long &delta) {
  if (DBG_LENSTR || DBG_GETMU_) {
    std::cout << "GETMU_: START" << std::endl;
  }

  //  return to_RR(2)*((to_RR(NumBits(delta) - 5) / to_RR(10)) + to_RR(6));
  //  return (to_RR(NumBits(delta) - 5) / to_RR(10)) + to_RR(6);
  RR mu = to_RR(6.85) + to_RR(10.62) * to_RR(NumBits(delta)) / to_RR(135);
  if (DBG_LENSTR || DBG_GETMU_) {
    std::cout << "GETMU_: FINISH" << std::endl;
  }
  return mu;
}

template <class U>
void RegulatorLenstraData<long, U>::bsgs_getentrysize(ZZ &entry_size,
                                                    bool nodist) {
  if (DBG_LENSTR || DBG_BSGSES) {
    std::cout << "BSGSES: START" << std::endl;
  }
  ZZ temp;

  entry_size = 3 * NTL_BITS_PER_LONG; // space for pointers
  if (DBG_BSGSES) {
    std::cout << "DBG_BSGSES: entry_size is " << entry_size << std::endl;
  }
  if (nodist) {
    if (DBG_BSGSES) {
      std::cout << "DBG_BSGSES: nodist is true" << std::endl;
    }
//     entry_size += SqrRoot(delta).size() * NTL_ZZ_NBITS; // a coeff
    entry_size += sizeof(delta) * NTL_ZZ_NBITS; // a coeff
    if (DBG_BSGSES) {
      std::cout << "DBG_BSGSES: entry_size is " << entry_size << std::endl;
    }
    entry_size += 8; // b coeff
    if (DBG_BSGSES) {
      std::cout << "DBG_BSGSES: entry_size is " << entry_size << std::endl;
    }
  } else {
    if (DBG_BSGSES) {
      std::cout << "DBG_BSGSES: nodist is false" << std::endl;
    }
    // temp = to_ZZ(1) << qo_distance<ZZ>::get_p();
    temp = to_ZZ(1) << 64;
    entry_size +=
//         temp.size() * NTL_ZZ_NBITS + 3 * SqrRoot(delta).size() * NTL_ZZ_NBITS;
        temp.size() * NTL_ZZ_NBITS + 3 * sizeof(delta) * NTL_ZZ_NBITS;
    if (DBG_BSGSES) {
      std::cout << "DBG_BSGSES: entry_size is " << entry_size << std::endl;
    }
  }

  // convert bits to bytes
  entry_size >>= 2;
  if (DBG_BSGSES) {
    std::cout << "DBG_BSGSES: End - entry_size is " << entry_size << std::endl;
  }
  if (DBG_LENSTR || DBG_BSGSES) {
    std::cout << "BSGSES: FINSH" << std::endl;
  }
}

// quadratic_order<ZZ>::init_prinlist
// Task: initializes the list of reduced principal ideals for BS-GS
// computations. Uses current contents if appropriate.

template <class U>
void RegulatorLenstraData<long, U>::init_prinlist(const ZZ &N, long l, U &s,
                                                long &M,
                                                QuadraticInfElement<long, U> &G) {

  if (DBG_LENSTR || DBG_IPLIST) {
    std::cout << "IPLIST: START" << std::endl;
  }

  if (DBG_IPLIST) {
    std::cout << "IPLIST: Begin" << std::endl;
    std::cout << "IPLIST: N is " << N << std::endl;
    std::cout << "IPLIST: l is " << l << std::endl;
  }

  long lsize, P;

  // compute B and s
  P = 3 + NumBits(SqrRoot(delta));
  s = to<U>((N + 2) * l);

  if (DBG_IPLIST) {
    std::cout << "IPLIST: s is " << s << std::endl;
    std::cout << "IPLIST: prin_list.no_of_elements() is " << prin_list.no_of_elements() << std::endl;
  }

  // initialize hash table
  if (prin_list.no_of_elements() > 0) {
    G.assign(prin_list.last_entry());
    if (DBG_IPLIST) {
      std::cout << "IPLIST: G.assign(prin_list.last_entry())" << std::endl;
      std::cout << "SHANKS: last_entry is " << prin_list.last_entry()
                << std::endl;
      std::cout << "SHANKS: G is " << G.get_qib() << " with distance "
                << G.get_distance() << std::endl;
    }
    l = prinlist_l;
    M = prinlist_M;

    if (prinlist_s > s) {
      s = prinlist_s;
    } else {
      prinlist_s = s;
    }
  } else {
    if (DBG_IPLIST) {
      std::cout << "IPLIST: before conv(lsize, N + P); " << std::endl;
      std::cout << "IPLIST: N is " << N << std::endl;
      std::cout << "IPLIST: P is " << P << std::endl;
      std::cout << "IPLIST: lsize is " << lsize << std::endl;
    }
    conv(lsize, N + P);
    lsize += 100;

    if (DBG_IPLIST) {
      std::cout << "IPLIST: before prin_list.initialize(lsize); " << std::endl;
      std::cout << "IPLIST: lsize is " << lsize << std::endl;
    }
    prin_list.initialize(lsize);

    if (DBG_IPLIST) {
      std::cout << "IPLIST: before G.assign_one(); " << std::endl;
    }
    G.assign_one();
    if (DBG_IPLIST) {
      std::cout << "IPLIST: G.assign_one()" << std::endl;
      std::cout << "SHANKS: G is " << G.get_qib() << " with distance "
                << G.get_distance() << std::endl;
    }
    prinlist_l = l;
    prinlist_s = s;
  }

  if (DBG_LENSTR || DBG_IPLIST) {
    std::cout << "IPLIST: FINISH" << std::endl;
  }
}

//
// quadratic_order<ZZ>::regulator_bsgs
// Task: Computes the regulator using an O(sqrt(R)) baby-step giant-step
// algorithm.
//

template <class U> void RegulatorLenstraData<long, U>::regulator_bsgs(ZZ &bound) {

  if (DBG_LENSTR || DBG_SHANKS) {
    std::cout << "SHANKS: START" << std::endl;
  }
  // initialize hash table
  ZZ K, B, N, entry_size;
  U u, s;
  RR mu;
  long l, M = 1;

  QuadraticInfElement<long, U> G{*quadratic_order};

  //   if(DBG_SHANKS) std::cout << "SHANKS: Begin" << std::endl;
  //   if(DBG_SHANKS) std::cout << "K is " << K << std::endl;
  //   if(DBG_SHANKS) std::cout << "delta is " << delta << std::endl;

  if (bound == 0) {
    // replacing get_bound for now...
    K = SqrRoot(delta);
    // get_bound(K);
  } else {
    K = bound;
  }

  if (DBG_SHANKS) {
    std::cout << "SHANKS: K is " << K << std::endl;
  }

  l = bsgs_getl(K, N, entry_size, mu, false);
  init_prinlist(N, l, s, M, G);

  B = N * l;

  // compute list of baby steps (distance < B)
  if (DBG_SHANKS) {
    std::cout << "SHANKS: G is " << G.get_qib() << " with distance "
              << G.get_distance() << std::endl;
    std::cout << "SHANKS: B is " << B << std::endl;
    std::cout << "SHANKS: Begin baby step list computation" << std::endl;
  }

  QuadraticInfElement<long, U> A{*quadratic_order}, C{*quadratic_order},
      AA{*quadratic_order};
  HashEntryReal<long, U> *F;

  A.assign_one();

  if (DBG_SHANKS) {
    std::cout << "SHANKS: Begin get_baby_steps" << std::endl;
    std::cout << "SHANKS: l == " << l << std::endl;
  }

  if (l == 1) {
    regulator = G.get_baby_steps(prin_list, B, A);
  } else {
    regulator = G.get_baby_steps(prin_list, B, A, l, M);
  }
  prinlist_M = M;

  if (DBG_SHANKS) {
    std::cout << "SHANKS: End get_baby_steps" << std::endl;
    std::cout << "SHANKS: prin_list is " << prin_list << std::endl;
    std::cout << "SHANKS: G is " << G.get_qib() << " with distance "
              << G.get_distance() << std::endl;
    std::cout << "SHANKS: regulator is " << regulator << std::endl;
  }

  if (!IsZero(regulator)) {
    set_case_type(" -- Found during bsgs baby-step -- ");
    Rbsgs = true;
  }

  if (IsZero(regulator)) {
    Rbsgs = false;
//     G.adjust(s);
    G.inverse_rho();
    G.inverse_rho();
    u = 2 * s;

    A = G;
     if (DBG_SHANKS){
       std::cout << "SHANKS: Begin giant_step" << std::endl;
       std::cout << "SHANKS: (" << G.get_qib().get_a() << ", "
                 << G.get_qib().get_b() << ", " << G.get_qib().get_c() << ") "
                 << G.get_distance() << std::endl;
     }
/*
     G.giant_step(G);

     if (DBG_SHANKS)
       std::cout << "SHANKS: (" << G.get_qib().get_a() << ", "
                 << G.get_qib().get_b() << ", " << G.get_qib().get_c() << ") "
                 << G.get_distance() << std::endl;
     if (DBG_SHANKS)
       std::cout << "SHANKS: End giant_step" << std::endl;
     sqr(G, G); makeshift square above
         if(DBG_SHANKS) std::cout << "SHANKS: Begin adjust, baby_step" <<
         std::endl;
     if (DBG_SHANKS)
       std::cout << "SHANKS: G distance is " << G.get_distance() << std::endl;
     G.adjust(u);
     if (DBG_SHANKS)
       std::cout << "SHANKS: G distance is " << G.get_distance() << std::endl;
     if (G.is_one())
       G.baby_step();
         if(DBG_SHANKS) std::cout << "SHANKS: End adjust, baby_step" <<
         std::endl;
     if (DBG_SHANKS)
       std::cout << "SHANKS: G distance is " << G.get_distance() << std::endl;*/
  }

  //   if(DBG_SHANKS) std::cout << "SHANKS: End baby step list computation" <<
  //   std::endl;
  //
  // compute giant steps until R is found or the bound is exceeded
  //
  //   if(DBG_SHANKS) std::cout << "SHANKS: Begin giant step traversal" <<
  //   std::endl;
  long i;

  //   if(DBG_SHANKS) std::cout << "SHANKS: s and u are " << s << " and " << u
  //   << std::endl;
  while (IsZero(regulator) && (bound == 0 || A.eval() <= bound)) {
    s += u;

    if (DBG_SHANKS) {
      std::cout << "SHANKS: Taking a giant_step" << std::endl;
      std::cout << "SHANKS: G is (" << G.get_qib().get_a() << ", "
                << G.get_qib().get_b() << ", " << G.get_qib().get_c() << ") "
                << G.get_distance() << std::endl;
    }

    if (DBG_SHANKS) {
      std::cout << "SHANKS: A is (" << A.get_qib().get_a() << ", "
                << A.get_qib().get_b() << ", " << A.get_qib().get_c() << ") "
                << A.get_distance() << std::endl;
    }
    A.giant_step(G);

    if (DBG_SHANKS) {
      std::cout << "SHANKS: A.giant_step(G) is (" << A.get_qib().get_a() << ", "
                << A.get_qib().get_b() << ", " << A.get_qib().get_c() << ") "
                << A.get_distance() << std::endl;
    }

    // mul(A, A, G); makeshift multiplication above
    // A.adjust(s);

    // search for A, rho_1(A), ..., rho_l(A) in the hash table
    if (DBG_SHANKS)
      std::cout << "SHANKS: A distance is " << A.get_distance() << std::endl;
    AA = A;
    if (DBG_SHANKS)
      std::cout << "SHANKS: AA distance is " << AA.get_distance() << std::endl;
    if (DBG_SHANKS)
      std::cout << "SHANKS: M  is " << M << std::endl;
    if (DBG_SHANKS)
      std::cout << "SHANKS: regulator is " << regulator << std::endl;
    for (i = 0; i < M && IsZero(regulator); ++i) {
      F = prin_list.search(AA.hash_real());
      if (F) {
        // found AA in the hash table!
        if (DBG_SHANKS) {
          std::cout << "SHANKS: found AA in the hash table! " << std::endl;
          std::cout << "SHANKS: AA distance is " << AA.get_distance()
                    << std::endl;
          std::cout << "SHANKS: F distance is " << F->get_d() << std::endl;
        }
        combine_BSGS(regulator, AA, F);
        if (DBG_SHANKS)
          std::cout << "SHANKS: regulator is " << regulator << std::endl;
        nuclose(C, FloorToZZ(regulator));
        C.adjust(regulator);
        regulator = C.get_distance();
        if (DBG_SHANKS)
          std::cout << "SHANKS: regulator is " << regulator << std::endl;
        set_case_type(" -- Found during bsgs giant-step -- ");
      }

//       else {
//         F = prin_list.search((AA.conjugate()).hash_real());
//         if (F) {
//           //found A^-1 in the hash table!
//           if (DBG_SHANKS) {
//             std::cout << "SHANKS: found AA^-1 in the hash table! " << std::endl;
//           }
//           combine_conj_BSGS(regulator, AA, F);
//           nuclose(C, FloorToZZ(regulator));
//           C.adjust(regulator);
//           regulator = C.get_distance();
//         } else if (i < M)
//           //AA.rho(); Swapping this out for AA.baby_step();
//           if (DBG_SHANKS)
//             std::cout << "SHANKS: AA distance is " << AA.get_distance()
//                       << std::endl;
//         AA.baby_step();
//         if (DBG_SHANKS)
//           std::cout << "SHANKS: AA distance is " << AA.get_distance()
//                     << std::endl;
//       }
    }
    if (DBG_SHANKS)
      std::cout << "SHANKS: regulator is " << regulator << std::endl;
  }
  if (DBG_SHANKS)
    std::cout << "SHANKS: End giant step traversal" << std::endl;
  if (!IsZero(regulator))
    Rconditional = false;

  if (DBG_SHANKS)
    std::cout << "SHANKS: End - regulator is " << regulator << std::endl;
  if (DBG_LENSTR || DBG_SHANKS) {
    std::cout << "SHANKS: FINISH" << std::endl;
  }
}

//
// quadratic_order<ZZ>::approximate_hR()
//
// Task:
//      returns an approximation of hR
//

template <class U> U RegulatorLenstraData<long, U>::approximate_hR() {

  if (DBG_LENSTR || DBG_APPRHR) {
    std::cout << "APPRHR: STARTING" << std::endl;
  }
  RR hR, FI;

  long n = get_optimal_Q_cnum();
  FI = l_function->approximateL1(n);

  if (DBG_APPRHR) {
    std::cout << "APPRHR: get_optimal_Q_cnum() returned " << n << std::endl;
    std::cout << "APPRHR: FI = l_function->approximateL1(n) is " << FI
              << std::endl;
  }

  if (quadratic_order->is_imaginary()) {
    // h = sqrt(delta) L / Pi
    hR = FI * SqrRoot(to_RR(-delta)) / ComputePi_RR();
    if (delta == -4)
      hR *= 2;
    if (delta == -3)
      hR *= 3;
  } else {
    // hR = sqrt(delta) L / 2
    hR = FI * SqrRoot(to_RR(delta)) / to_RR(2);
  }

  if (DBG_LENSTR || DBG_APPRHR) {
    std::cout << "APPRHR: FINISHED" << std::endl;
  }
  return to<U>(hR);
}

template <class U>
void RegulatorLenstraData<long, U>::optimize_K(ZZ &bound, const ZZ &S,
                                             const ZZ &N, long l) {
  if (DBG_LENSTR || DBG_OPTIMK) {
    std::cout << "OPTIMK: STARTING" << std::endl;
  }
  RR K = to_RR(bound);
  RR mu = get_mu(delta);
  RR rS = to_RR(S);
  RR rN = to_RR(N);
  RR rn;

  if (DBG_OPTIMK) {
    std::cout << "K0 = " << bound << endl;
    std::cout << "mu = " << mu << endl;
    std::cout << "S = " << S << endl;
    std::cout << "N = " << N << endl;
    std::cout << "l = " << l << endl;
    std::cout << endl;
  }

  if (parallel)
    rn = to_RR(parallel);
  else
    set(rn);

  RR F, dF;

  F = func(K, rS, rN, mu, rn, l);
  dF = dfunc(K, rS, rN, mu, rn, l);
  if (DBG_OPTIMK) {
    std::cout << "K = " << K << ", F(K) = " << F << ", F'(K) = " << dF
              << std::endl;
  }

  while (abs(F) > to_RR(0.00001) && K - F / dF > 0) {
    K = K - F / dF;
    F = func(K, rS, rN, mu, rn, l);
    dF = dfunc(K, rS, rN, mu, rn, l);
  }
  bound = FloorToZZ(K);
  if (DBG_LENSTR || DBG_OPTIMK) {
    std::cout << "OPTIMK: FINISHED" << std::endl;
  }
}

template <class U>
void RegulatorLenstraData<long, U>::combine_BSGS(
    U &dist, const QuadraticInfElement<long, U> &DD,
    const HashEntryReal<long, U> *F) {

  dist = DD.get_distance() - F->get_d();
}

template <class U>
void RegulatorLenstraData<long, U>::combine_conj_BSGS(
    U &dist, const QuadraticInfElement<long, U> &DD,
    const HashEntryReal<long, U> *F) {

  dist = F->get_d() - DD.get_distance() - log(to<U>(DD.get_qib().get_a()));
}

//
// lower_bound_hR()
// Task: returns a lower bound of hR such that L < hR < 2L
//
template <class U> RR RegulatorLenstraData<long, U>::lower_bound_hR() {

  if (DBG_LENSTR || DBG_LOBOHR) {
    std::cout << "LOBOHR: STARTING" << std::endl;
  }
  RR hR;
  RR temp, FI;

  // START: Temporary variables needed (previously declared in ANTL-Import's
  // quadratic_order
  bool unconditional = false;
  int info = 0;
  bool use_tables = false;
  // FINISH: Temporary variables needed (previously declared in ANTL-Import's
  // quadratic_order

  if (unconditional) {
    if (info > 1) {
      std::cout << "Lower bound hR with L(0,X)" << std::endl;
    }

    // compute hR apporximation

    if (quadratic_order->is_imaginary()) {
      if (use_tables)
        FI = l_function->approximateL0_ImaginaryNumberField_table(
            log(sqrt(double(2))));
      else
        FI = l_function->approximateL0_ImaginaryNumberField(
            log(sqrt(double(2))));

      if (info > 1) {
        std::cout << "L(0,X) approx " << FI << std::endl;
      }

      temp = FI;
      if (delta == -4)
        temp *= 2;
      if (delta == -3)
        temp *= 3;

    } else {
      if (use_tables)
        FI = l_function->approximateL0_RealNumberField_table(
            log(sqrt(double(2))));
      else
        FI = l_function->approximateL0_RealNumberField(log(sqrt(double(2))));

      if (info > 1) {
        std::cout << "L(0,X) approx " << FI << std::endl;
      }
      temp = FI;
    }
  } else {
    if (info > 1) {
      std::cout << "Lower bound hR with L(1,X)" << std::endl;
    }

    if (use_tables)
      FI = l_function->approximateL1_table();
    else {
      long n = get_optimal_Q();
      if (info > 2)
        std::cout << "using n = " << n << std::endl;
      FI = l_function->approximateL1(n);
    }

    if (info > 1) {
      std::cout << "L(1,X) approx " << FI << std::endl;
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

  hR = temp / SqrRoot(to_RR(2));

  if (info > 1)
    std::cout << "Lower bound = " << hR << std::endl;

  if (DBG_LENSTR || DBG_LOBOHR) {
    std::cout << "LOBOHR: FINISH" << std::endl;
  }
  return hR;
}

//
// quadratic_order<ZZ>::get_optimal_Q()
//
// Task: returns a value of Q from a pre-computed table which will compute h*
// such that h* < h < 2h*.
//

template <class U> long RegulatorLenstraData<long, U>::get_optimal_Q() {
  if (DBG_LENSTR || DBG_GETOPQ) {
    std::cout << "GETOPQ: STARTING" << std::endl;
  }
  long Dlog;
  ZZ temp;

  temp = FloorToZZ(log10(to_RR(abs(delta))));
  conv(Dlog, temp);

  if ((Dlog / 5) > 19)
    return generate_optimal_Q();

  long optimal_Q = OQvals[Dlog / 5];

  if (DBG_LENSTR || DBG_GETOPQ) {
    std::cout << "GETOPQ: FINISH" << std::endl;
  }
  return optimal_Q;
}

//
// quadratic_order<ZZ>::generate_optimal_Q()
// Task: computes the optimal value of Q for computing h* such that h* < h < 2h*
//

template <class U> long RegulatorLenstraData<long, U>::generate_optimal_Q() {
  if (DBG_LENSTR || DBG_GENOPQ) {
    std::cout << "GENOPQ: STARTING" << std::endl;
  }
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

  if (DBG_LENSTR || DBG_GENOPQ) {
    std::cout << "GENOPQ: FINISHED" << std::endl;
  }
  return OQ;
}

//
// quadratic_order<T>::find_hstar
//
// Task:
//      finds the regulator unconditionally given a multiple.  If S is the
//      multiply, requires S^1/3 operations.
//

template <class U>
void RegulatorLenstraData<long, U>::find_hstar(ZZ &hstar, const U &S,
                                             const ZZ &Pmax) {
  if (DBG_LENSTR || DBG_FHSTAR) {
    std::cout << "FHSTAR: STARTING" << std::endl;
  }

  set(hstar);

  if (DBG_FHSTAR) {
    std::cout << "FHSTAR: hstar = " << hstar << std::endl;
    std::cout << "FHSTAR: to_RR(S) = " << to_RR(S) << std::endl;
    std::cout << "FHSTAR: S = " << S << std::endl;
    std::cout << "FHSTAR: Pmax = " << Pmax << std::endl;
  }

  // Generate a list of primes
  PrimeSeq prime_seq;
  long prime = prime_seq.next();

  // Preparation for
  QuadraticInfElement<long, U> target_qie{*quadratic_order};
  U target_distance;

  // Compute the power of each prime in the factorization of hstar
  if (DBG_FHSTAR) {
    std::cout << "FHSTAR: Testing each prime" << std::endl;
  }
  while (prime <= Pmax) {
    if (DBG_FHSTAR) {
      std::cout << "FHSTAR: prime is " << prime << std::endl;
    }
    int power = 1;
    target_distance = S / to<U>(prime);
    nuclose(target_qie, FloorToZZ(target_distance));
    target_qie.adjust(target_distance);

    if (DBG_FHSTAR) {
      std::cout << "FHSTAR: target_distance was " << target_distance
                << std::endl;
      std::cout << "FHSTAR: target_qie was " << target_qie.get_qib()
                << std::endl;
      std::cout << "FHSTAR: target_qie distance was "
                << target_qie.get_distance() << std::endl;
    }

    // Increase prime power until corresponding distance is no longer a multiple
    // of the regulator
    if (DBG_FHSTAR) {
      std::cout << "FHSTAR: Finding power of current prime" << std::endl;
    }
    while (target_qie.is_one() &&
           abs(target_qie.get_distance() - target_distance) < 0.01) {
      if (DBG_FHSTAR) {
        std::cout << "FHSTAR: power is " << power << std::endl;
      }
      power++;
      target_distance = target_distance / to<U>(prime);
      nuclose(target_qie, FloorToZZ(target_distance));
      target_qie.adjust(target_distance);

      if (DBG_FHSTAR) {
        std::cout << "FHSTAR: target_distance was " << target_distance
                  << std::endl;
        std::cout << "FHSTAR: target_qie was " << target_qie.get_qib()
                  << std::endl;
        std::cout << "FHSTAR: target_qie distance was "
                  << target_qie.get_distance() << std::endl;
      }
    }
    if (DBG_FHSTAR) {
      std::cout << "FHSTAR: hstar was " << hstar << std::endl;
      std::cout << "FHSTAR: prime was " << prime << std::endl;
      std::cout << "FHSTAR: power was " << power - 1 << std::endl;
    }
    hstar *= FloorToZZ(pow(double(prime), double(power - 1)));
    if (DBG_FHSTAR) {
      std::cout << "FHSTAR: hstar is " << hstar << std::endl;
    }
    prime = prime_seq.next();
  }

  // Using hstar, computer the regulator
  target_distance = S / to<U>(hstar);
  nuclose(target_qie, FloorToZZ(target_distance));
  target_qie.adjust(target_distance);
  regulator = target_qie.get_distance();

  if (DBG_FHSTAR) {
    std::cout << "FHSTAR: hstar = " << hstar << std::endl;
    std::cout << "FHSTAR: R = " << regulator << std::endl;
  }
  if (DBG_LENSTR || DBG_FHSTAR) {
    std::cout << "FHSTAR: FINISHED" << std::endl;
  }
}

#endif
