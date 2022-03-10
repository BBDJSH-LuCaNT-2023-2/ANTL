#ifndef REGULATOR_LENSTRA_COMPUTE_H
#define REGULATOR_LENSTRA_COMPUTE_H

#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstraData.hpp>
#include <ANTL/common.hpp>

NTL_CLIENT
using namespace ANTL;

namespace ANTL {

template <class T> class RegulatorLenstraCompute {
public:
  RR regulator_lenstra(RegulatorLenstraData<T> &RegulatorLenstraData);

private:
  ZZ estimate_hR_error(RegulatorLenstraData<T> &RegulatorLenstraData);
};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraCompute_impl.hpp file.

// quadratic_order<T>::regulator_shanks
// Task: Computes the regulator using an O(D^1/5) baby-step giant-step algorithm
// of Shanks.

template <class T>
RR RegulatorLenstraCompute<T>::regulator_lenstra(RegulatorLenstraData<T> &RegulatorLenstraData) {

  //
  // initialize hash table
  //

  RR regulator;
  ZZ K, N, B, entry_size, u, s, s2;
  long l = 1, M = 1;
  RR mu;
  qi_pair<T> AA, G;

  qo_distance<T> S;
  clear(S);

  K = estimate_hR_error(RegulatorLenstraData.get_quadratic_order()) >> 1;
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

  qi_pair<T> A, C, CC, D, DD, GG;
  qo_hash_entry_real<T> *F;

  if (IsZero(S)) {

    A.assign_one();
    if (l == 1)
      regulator = G.get_baby_steps(prin_list, B, A);
    else
      regulator = G.get_baby_steps(prin_list, B, A, l, M);

    if (!IsZero(R)) {
      nuclose(C, R.eval());
      regulator = C.get_distance();
      Rbsgs = true;
    }

    prinlist_M = M;

    if (IsZero(regulator)) {
      G.adjust(s);

      u = (s << 1);
      nudupl(G, G);
      G.adjust(u);
      if (G.is_one())
        G.rho();

      GG = conjugate(G);
      while (abs(GG.get_distance().eval()) <= u)
        GG.inverse_rho();
      while (abs(GG.get_distance().eval()) > u)
        GG.rho();

      s2 = s = AA.get_distance().eval();
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
        F = prin_list.search((-CC).hash_real());
        if (F) {
          // found CC^-1 in the hash table!

          combine_conj_BSGS(S, CC, F);
        } else {
          F = prin_list.search(DD.hash_real());
          if (F) {
            // found DD in the hash table!

            combine_BSGS(S, DD, F);
          } else {
            F = prin_list.search((-DD).hash_real());
            if (F) {
              // found DD^-1 in the hash table!

              combine_conj_BSGS(S, DD, F);
            } else if (i < M) {
              CC.rho();
              DD.rho();
            }
          }
        }
      }
    }

    if (IsZero(S)) {
      s += u;
      nucomp(C, C, G);
      C.adjust(s);

      s2 -= u;
      nucomp(D, D, GG);
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
    RR temp = log(to_RR(S.eval())) * to_RR(2) / to_RR(3);

    conv(B, ceil(exp(temp)));

    l = bsgs_getl(B, N, entry_size, mu, true);
    optimize_K(B, S.eval(), N, l);

    regulator_bsgs(B);

    ZZ hstar, Pmax;

    if (info > 1 && !IsZero(regulator)) {
      hstar = S.eval() / R.eval();

      cout << "Found R = " << regulator << endl;
      cout << "h* = " << hstar << endl;
    } else {
      Pmax = 1 + S.eval() / B;
      find_hstar(hstar, S.eval(), Pmax, t1);
    }

    nuclose(C, R.eval());
    regulator = C.get_distance();
  }
}

// quadratic_order<ZZ>::estimate_hR_error(long n)
// Task: returns L such that |hR - hR'| < exp(L)^2

template <> ZZ RegulatorLenstraCompute<ZZ>::estimate_hR_error(RegulatorLenstraData<ZZ> &RegulatorLenstraData) {
  ZZ err;
  RR Aval, Fval, temp;

  if (RegulatorLenstraData.get_l_function().terms_used(1) == 0)
    return ZZ::zero();

  long n = get_optimal_Q_cnum();
  RR FI = RegulatorLenstraData.get_l_function().approximateL1(n);

  Aval = RegulatorLenstraData.get_l_function().calculate_L1_error(Delta, l_function.terms_used(1));
  Fval = exp(Aval) - 1;
  temp = 1 - exp(-Aval);
  if (temp > Fval)
    Fval = temp;

  if (RegulatorLenstraData.get_quadratic_order().is_imaginary()) {
    // h = sqrt(Delta) L / Pi
    Fval *= FI * SqrRoot(to_RR(-Delta)) / ComputePi_RR();
    if (Delta == -4)
      temp *= 2;
    if (Delta == -3)
      temp *= 3;
  } else {
    // hR = sqrt(Delta) L / 2
    Fval *= FI * SqrRoot(to_RR(Delta)) / 2;
  }

  err = CeilToZZ(Fval * log(to_RR(2))) >> 1;
  return err;
}
} // namespace ANTL
#endif
