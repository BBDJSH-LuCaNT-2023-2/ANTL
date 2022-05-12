#ifndef CLASSGROUP_BSGS_H
#define CLASSGROUP_BSGS_H

#include <ANTL/common.hpp>
#include <unordered_set>

#include <ANTL/HashTable/HashEntryInt.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>

#include <ANTL/LinearAlgebra/Smith.hpp>

#include <ANTL/Quadratic/QuadraticClassGroupElement.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

template<>
struct std::hash<QuadraticIdealBase<ZZ>>
{
    std::size_t operator()(QuadraticIdealBase<ZZ> const& qib) const noexcept
    {
        std::size_t h1 = std::hash<int>{}(to<int>(qib.get_a()));
        std::size_t h2 = std::hash<int>{}(to<int>(qib.get_b()));
        return h1 ^ (h2 << 1); // or use boost::hash_combine
    }
};

NTL_CLIENT;
using namespace ANTL;

namespace ANTL {

// DBG_CONSTANTS
bool DBG_CGBGRL = true;

void decode_vector(mat_ZZ &Bmat, const ZZ &Bjj, const ZZ &r, const ZZ &q,
                   vec_ZZ &Rvec, vec_ZZ &Qvec, long nR, long nQ);

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class T> class ClassGroupBSGSReal {
private:
  QuadraticOrder<T> *quadratic_order;

  double regulator;

  T delta;

  // factor base for computing CL
  QuadraticClassGroupElement<T> *fact_base;

  // indice of fact base elements that
  long *contributors;

  // number of primes used for BSGS
  long num_prime_ideals;

  // class number
  ZZ h;

  // invariants of CL
  vec_ZZ CL;

  // size of fact_base
  long numFB;

  // transformation matrix
  mat_ZZ U_mat;

public:
  ClassGroupBSGSReal(QuadraticOrder<T> *quadratic_order);

  void cg_bsgs_real(const ZZ &hstar);

  void set_regulator(double &ext_regulator){regulator = ext_regulator;}

  vec_ZZ get_class_group(){return CL;}

private:
  T get_dist_mod(const T &Delta) { return CeilToZZ(log(to_RR(Delta))); }

  void get_next_prime(QuadraticClassGroupElement<T> &G);

  bool is_principal(const QuadraticClassGroupElement<T> &G);
};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class T>
ClassGroupBSGSReal<T>::ClassGroupBSGSReal(QuadraticOrder<T> *quadratic_order_arg) {

  quadratic_order = quadratic_order_arg;
  delta = quadratic_order->get_discriminant();

//   parallel = false;
//   Rbsgs = false;        // true if R was computed using BSGS
//   Rconditional = false; // true if correctness of R relies on ERH
}

template <class T>
void ClassGroupBSGSReal<T>::cg_bsgs_real(const ZZ &hstar) {
  if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 1 - Initialize variables" << std::endl;}

  ZZ y, usqr, det, Bjj, temp, curr_index;
  ZZ s, u, r, q, Bj, hS;
  long i, j, k, upper, crank, numRpr, numQ, idx;

  QuadraticClassGroupElement<T> G{*quadratic_order}, A{*quadratic_order},
      B{*quadratic_order}, C{*quadratic_order}, D{*quadratic_order},
      E{*quadratic_order}, HI{*quadratic_order}, Gq{*quadratic_order},
      GBj{*quadratic_order}, RHO{*quadratic_order};

  mat_ZZ Bmat, junk;
  vec_ZZ Rvec, Qvec;

  vec_long Rideals;
  long numRideals = 0;
  QuadraticInfElement<T, double> F{*quadratic_order}, RHOdist{*quadratic_order};
  ZZ sqReg, Fstep, curr_dist, diff;
  bool done;

  HashEntryInt<T, ZZ> *Inode;
  IndexedHashTable<HashEntryInt<T, ZZ>> QT, RT;
  //  RR nFI = to<RR>(hstar) * SqrRoot (to_RR (2));

  if (DBG_CGBGRL)
    std::cout << "CG_BSGS_REAL:  hstar = " << hstar << std::endl;

  if (hstar > 1) {
    Rvec.kill();
    Qvec.kill();
    Rideals.SetLength(1024);

    // get giant step for equivalence testing
    SqrRoot(sqReg, regulator);
    diff = get_dist_mod(delta);
    F.assign_one();

    if (-0.000000001 < regulator - to<double>(sqReg + get_dist_mod(delta)))
      nuclose(F, sqReg);

    if (F.get_qib().IsOne())
      Fstep = regulator;
    else
      Fstep = F.get_distance();

    if (DBG_CGBGRL) {std::cout << "sqReg = " << sqReg << std::endl;}
    if (DBG_CGBGRL) {std::cout << "regulator = " << regulator << std::endl;}
    if (DBG_CGBGRL) {std::cout << "mod = " << get_dist_mod(delta) << std::endl;}
    if (DBG_CGBGRL) {std::cout << "Fstep = " << Fstep << std::endl;}
//     if (DBG_CGBGRL)
//       std::cout << "F = " << F << std::endl;

    if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 2 - Initialize sets R and Q" << std::endl;}
    // initialize sets R and Q
    temp = SqrRoot(hstar << 1);
    upper = to<long>(temp);

    QT.initialize(upper << 1);
    B.assign_one();
    QT.hash(B.hash_int(to<ZZ>(0)));

    RT.initialize(upper << 1);
    //    RT.hash (B.hash_int (to<S>(0)));

    Rideals[0] = 0;
    numRideals = 1;

    RHO.assign(B);
    RHOdist.assign(B);
    do {
      RT.hash(RHO.hash_int(to<ZZ>(0)));
      RHOdist.baby_step();
      RHO.assign(RHOdist.get_qib());
      if (F.get_qib().IsOne())
        done = (RHO.IsOne());
      else
        done = (RHOdist.get_distance() > sqReg + diff || RHO.IsOne());
    } while (!done);
  }

  if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 3 - Initialize fact_base and contributors" << std::endl;}
  fact_base = new QuadraticClassGroupElement<T>[100];
  contributors = new long[100];

  num_prime_ideals = 0;
  crank = 0;
  NTL::set(det);
  NTL::set(hS);
  G.assign_one();

  if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4 - (hS < hstar) loop" << std::endl;}
  // while h < hstar, we only have a subgroup
  while (hS < hstar) {
    // get prime ideal
    if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.1 - the do while(is_principal(G)) loop" << std::endl;}
    do {
      get_next_prime(G);
      ++num_prime_ideals;
    } while (is_principal(G));

    fact_base[crank] = G;
    contributors[crank] = crank;

    if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.2 - set initial step-width" << std::endl;}
    // set initial step-width
    u = y = 1;
    s = 1;

    A.assign_one();
    nupower_real(C, G, to<ZZ>(u));
    B.assign(C);
    HI.assign(conjugate(G));

    if (DBG_CGBGRL)
      std::cout << "\nG = " << G << std::endl;
    if (DBG_CGBGRL)
      std::cout << "G^-1 = " << HI << std::endl;
    if (DBG_CGBGRL)
      std::cout << "u = " << u << ", y = " << y << std::endl;
    if (DBG_CGBGRL)
      std::cout << "B = " << B << std::endl;
    if (DBG_CGBGRL)
      std::cout << "C = " << C << std::endl;
    if (DBG_CGBGRL)
      std::cout << "crank = " << crank << std::endl;
    if (DBG_CGBGRL)
      std::cout << "Rvec = " << Rvec << std::endl;
    if (DBG_CGBGRL)
      std::cout << "Qvec = " << Qvec << std::endl;

    curr_index = numRideals;
    numQ = QT.no_of_elements();

    ::clear(Bjj);

    if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.3 - check whether current ideal is in previously generated subgroup" << std::endl;}
    // check whether current ideal is in previously generated subgroup
    for (i = 0; i < numQ && ::IsZero(Bjj); ++i) {

      if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.3.0 - assign;" << std::endl;}
      D.assign(QT[i]);
      if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.3.1 - mul(E, D, G);" << std::endl;}
//       nucomp_real(E, D, G);
      mul(E, D, G);

      if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.3.2 - test for all RHO equivalent to E" << std::endl;}
      // test for all RHO equivalent to E
      if (F.get_qib().IsOne()) {
        // entire equivalence classes are stored in RT
        Inode = RT.search(E.hash_int(to<ZZ>(0)));
        if (Inode) {
          NTL::set(Bjj);
        }
      } else {
        // use baby-step giant-step for equivalence testing
        // baby steps are all in RT
        RHOdist.assign(E);
        ::clear(curr_dist);
        do {
          RHO.assign(RHOdist.get_qib());
          Inode = RT.search(RHO.hash_int(to<ZZ>(0)));
          if (Inode) {
            NTL::set(Bjj);
          }

          if (::IsZero(Bjj)) {
            curr_dist += Fstep;
//             nucomp(RHOdist, RHOdist, F);
            RHOdist.giant_step(F);
            RHOdist.adjust(curr_dist);
          }
        } while (RHOdist.get_distance() <= regulator && ::IsZero(Bjj));
      }
    }

    if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 4.3 - while (::IsZero(Bjj))" << std::endl;}
    while (::IsZero(Bjj)) {

      if (DBG_CGBGRL)
        std::cout << "\nBABY" << std::endl;
      if (DBG_CGBGRL)
        std::cout << "s = " << s << ", u = " << u << ", G = " << G
                  << ", numRpr = " << RT.no_of_elements()
                  << ", numRideals = " << numRideals
                  << ", curr_index = " << curr_index << std::endl;

      // compute more baby steps
      for (r = s; r <= u && ::IsZero(Bjj); ++r) {
//         nucomp_real(A, A, HI);
        mul(A, A, HI);

        if (DBG_CGBGRL)
          std::cout << "\nr = " << r << ", A = " << A << std::endl;

        for (k = 0; k < numRideals; ++k) {
          D.assign(RT[Rideals[k]]);
//           nucomp_real(E, D, A);
          mul(E, D, A);

          if (curr_index == Rideals.MaxLength())
            Rideals.SetLength(Rideals.MaxLength() << 1);

          Rideals[to<long>(curr_index)] = RT.no_of_elements();

          if (DBG_CGBGRL)
            std::cout << "BABY (NS):  r = " << r << ", E = " << E << std::endl;

          RHO.assign(E);
          RHOdist.assign(E);
          do {
            RT.hash(RHO.hash_int(curr_index));
            RHOdist.baby_step();
            RHO.assign(RHOdist.get_qib());
            if (F.get_qib().IsOne())
              done = (RHO == E);
            else
              done = (RHOdist.get_distance() > sqReg + diff || RHO == E);
          } while (!done);

          ++curr_index;
        }
      }

      if (DBG_CGBGRL)
        std::cout << "Got baby"
                  << ", numRT = " << RT.no_of_elements()
                  << ", curr_index = " << curr_index << ", numQ = " << numQ
                  << std::endl;
      if (DBG_CGBGRL)
        std::cout << "\nu = " << u << ", u^2 = " << u * u << std::endl;

      // compute giant steps to u^2
      usqr = u * u;
      while ((y < usqr) && ::IsZero(Bjj)) {

        if (DBG_CGBGRL)
          std::cout << "\nGIANT:  y = " << y << ", u^2 = " << usqr
                    << ", B = " << B << std::endl;

        // search
        for (i = 0; i < numQ && ::IsZero(Bjj); ++i) {
          E.assign(QT[i]);
//           nucomp_real(D, E, B);
          mul(D, E, B);

          // test for all RHO equivalent to E
          if (F.get_qib().IsOne()) {
            // entire equivalence classes are stored in RT
            Inode = RT.search(D.hash_int(to<ZZ>(0)));
            if (Inode && !D.IsOne()) {

              if (DBG_CGBGRL)
                std::cout << "FOUND GIANT" << std::endl;

              r = Inode->get_d();
              q = i;
              Bjj = y + r / numRideals;

              if (DBG_CGBGRL)
                std::cout << "r = " << r << ", q = " << q << ", Bjj = " << Bjj
                          << std::endl;

              if (Bjj > 1) {
                r %= to<ZZ>(numRideals);
                decode_vector(Bmat, to<ZZ>(Bjj), to<ZZ>(r), to<ZZ>(q), Rvec,
                              Qvec, numRideals, numQ);

                // PART_OF_OLD_DEBUG: test_relation_cg(fact_base, crank, Bmat);
              }
            }
          } else {
            // use baby-step giant-step for equivalence testing
            // baby steps are all in RT
            RHOdist.assign(D);
            ::clear(curr_dist);
            do {
              RHO.assign(RHOdist.get_qib());
              Inode = RT.search(RHO.hash_int(to<ZZ>(0)));
              if (Inode && !D.IsOne()) {

                if (DBG_CGBGRL)
                  std::cout << "FOUND GIANT" << std::endl;

                r = Inode->get_d();
                q = i;
                Bjj = y + r / numRideals;

                if (DBG_CGBGRL)
                  std::cout << "r = " << r << ", q = " << q << ", Bjj = " << Bjj
                            << std::endl;

                if (Bjj > 1) {
                  r %= to<ZZ>(numRideals);
                  decode_vector(Bmat, to<ZZ>(Bjj), to<ZZ>(r), to<ZZ>(q), Rvec,
                                Qvec, numRideals, numQ);
                  // PART_OF_OLD_DEBUG: test_relation_cg(fact_base, crank,
                  // Bmat);
                }
              }

              if (::IsZero(Bjj)) {
                curr_dist += Fstep;
//                 nucomp(RHOdist, RHOdist, F);
                RHOdist.giant_step(F);
                RHOdist.adjust(curr_dist);
              }
            } while (RHOdist.get_distance() <= regulator && ::IsZero(Bjj));
          }
        }

        if (::IsZero(Bjj)) {
          // not found, take another giant step
          y += u;
//           nucomp_real(B, B, C);
          mul(B, B, C);
        }
      }

      if (DBG_CGBGRL) {std::cout << "Got giant" << std::endl;}

      // double u
      if (::IsZero(Bjj)) {
        s = u + 1;
        u <<= 1;
//         nudupl_real(C, C);
        mul(C, C, C);
      }
      if (DBG_CGBGRL) {std::cout << "Just after Got giant" << std::endl;}
    }

    if (!::IsOne(Bjj)) {

      hS *= Bjj;
      //	  nFI /= to<RR> (Bjj);
      ++crank;
    }

    if (DBG_CGBGRL)
      std::cout << "h = " << hS << std::endl;

    if (hS < hstar) {
      Bj = to<ZZ>(CeilToZZ(sqrt(to<RR>(Bjj))));

      // compute new R' (remove entries with too large exponents)
      det *= Bj;
      idx = to<long>(det);
      numRpr = RT.no_of_elements();

      if (DBG_CGBGRL)
        std::cout << "REMOVING..." << std::endl;

      i = numRpr - 1;
      while (RT[i].get_d() >= idx) {
        RT.remove_from(i);
        --i;
      }
      numRideals = idx;

      if (DBG_CGBGRL)
        std::cout << "DONE" << std::endl;

      if (!::IsOne(Bjj)) {
        Rvec.SetLength(crank);
        Qvec.SetLength(crank);
        Rvec[crank - 1] = Bj;
        Qvec[crank - 1] = 1;

        if (DBG_CGBGRL)
          std::cout << "SET LEN" << std::endl;

        // compute new Q
        numQ = QT.no_of_elements();
        curr_index = numQ;
        nupower_real(GBj, G, to<ZZ>(Bj));
        Gq.assign(GBj);

        if (DBG_CGBGRL)
          std::cout << "NEW Q:  Bj = " << Bj << std::endl;

        for (i = 1; i < Bj; ++i) {
          if (Bjj <= (i * Bj))
            break;
          for (k = 0; k < numQ; ++k) {
            E.assign(QT[k]);
//             nucomp_real(D, E, Gq);
            mul(D, E, Gq);
            QT.hash(D.hash_int(curr_index));
            ++curr_index;
          }
//           nucomp_real(Gq, Gq, GBj);
          mul(Gq, Gq, GBj);
          ++Qvec[crank - 1];
        }
        if (DBG_CGBGRL)
          std::cout << "DONE Q" << std::endl;
      }
    }
  }

  if (DBG_CGBGRL) {std::cout << "CGBGRL: STEP 5 - compute structure" << std::endl;}
  // compute structure
  if (DBG_CGBGRL)
    std::cout << "COMPUTING STRUCTURE" << std::endl;

  h = to<ZZ>(hS);
  CL.SetLength(1);
  if (h == 1) {
    CL[0] = 1;
    fact_base[0].assign_one();
    contributors[0] = 0;
    numFB = 1;
  } else {
    // delete all rows and columns with diagonal 1
    numFB = crank;

    mat_ZZ SNF;
    SNF = SmithNormalForm(Bmat, U_mat, junk);

    i = 0;
    for (j = 0; j < crank; ++j) {
      Bjj = to<ZZ>(SNF[j][j]);
      if (!::IsOne(Bjj)) {
        CL.SetLength(i + 1);
        CL[i++] = to<ZZ>(Bjj);
      }
    }
  }
}

//
// quadratic_order<T>::get_next_prime
// Task: Computes the next smallest prime ideal (in terms of norm)
//
template <class T> void ClassGroupBSGSReal<T>::get_next_prime(QuadraticClassGroupElement<T> &G) {
  static ZZ X, q;
  T P;

  if (G.IsOne()) {
    q = CARDINALITY<T>();
    X = q - 1;
  }

  do {
    do {
      ++X;
      get_poly_modq(P, X, q);
    } while (!::IsOne(LeadCoeff(P)) || !DetIrredTest(P));
  } while (!G.assign_prime(P));
}

template <>
void ClassGroupBSGSReal<ZZ>::get_next_prime(QuadraticClassGroupElement<ZZ> &G) {
  static PrimeSeq PS;

  if (G.IsOne())
    PS.reset(0);

  long p = PS.next();

  while (!G.assign_prime(to_ZZ(p)))
    p = PS.next();
}

template <class T>
bool ClassGroupBSGSReal<T>::is_principal(
    const QuadraticClassGroupElement<T> &G) {
  double sqrt_regulator = sqrt(regulator);

  QuadraticInfElement<T, double> baby_step_list_generator{*quadratic_order};

  std::unordered_set<QuadraticIdealBase<T>> baby_step_list = {};
  baby_step_list_generator.baby_step();

  while (baby_step_list_generator.get_distance() < sqrt_regulator &&
         baby_step_list.count(baby_step_list_generator.get_qib()) == 0) {
    baby_step_list.insert(baby_step_list_generator.get_qib());
    baby_step_list_generator.baby_step();
  }

  if (baby_step_list.count(G) != 0) {
    return true;
  }

  else {

    QuadraticClassGroupElement<T> giant_step{*quadratic_order}, G_copy{*quadratic_order};
    G_copy.assign(G);
    giant_step.assign(G);


    for(int i = 0; i <= CeilToZZ(to_RR(sqrt_regulator)); i++){
      mul(G_copy, G_copy, giant_step);
      if(baby_step_list.count(G) != 0){
        return true;
      }
    }
  }
  return false;
}



void decode_vector(mat_ZZ &Bmat, const ZZ &Bjj, const ZZ &r, const ZZ &q,
                   vec_ZZ &Rvec, vec_ZZ &Qvec, long nR, long nQ) {
  long i, j, rrank;
  ZZ temp, nRR, nQQ, rr, qq;
  vec_ZZ Bj;
  mat_ZZ tmat;

  rr = r;
  qq = q;
  nRR = nR;
  nQQ = nQ;

  // compute powers in index vector
  rrank = Bmat.NumCols();

  if (rrank == 0) {
    Bmat.SetDims(1, 1);
    Bmat[0][0] = Bjj;
  } else {
    Bj.SetLength(rrank + 1);

    Bj[rrank] = Bjj;
    for (i = rrank - 1; i >= 0; --i) {
      nRR /= Rvec[i];
      nQQ /= Qvec[i];
      Bj[i] = (rr / nRR) + (qq / nQQ) * Rvec[i];
      rr %= nRR;
      qq %= nQQ;
    }

    // new row and column
    ++rrank;
    tmat.SetDims(rrank, rrank);
    for (i = 0; i < rrank - 1; ++i)
      for (j = 0; j < rrank - 1; ++j)
        tmat[i][j] = Bmat[i][j];

    for (i = 0; i < rrank; ++i)
      tmat[i][rrank - 1] = Bj[i];

    Bmat.kill();
    Bmat.SetDims(rrank, rrank);
    Bmat = tmat;
  }
}

} // namespace ANTL
#endif
