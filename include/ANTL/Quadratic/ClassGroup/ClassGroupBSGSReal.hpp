#ifndef CLASSGROUP_BSGS_H
#define CLASSGROUP_BSGS_H

#include <ANTL/common.hpp>

#include <ANTL/HashTable/HashEntryInt.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>

#include <ANTL/Quadratic/QuadraticClassGroupElement.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

NTL_CLIENT;
using namespace ANTL;

namespace ANTL {

// DBG_CONSTANTS
bool DBG_CGBGRL = false;

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class T, class U> class ClassGroupBSGSReal {
private:
  U regulator;

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
  void cg_bsgs_real(const U &hstar);

private:
};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class T, class U>
void ClassGroupBSGSReal<T, U>::cg_bsgs_real(const U &hstar) {
  U y, usqr, det, Bjj, temp, curr_index;
  U s, u, r, q, Bj, hS;
  long i, j, k, upper, crank, numRpr, numQ, idx;

  QuadraticClassGroupElement<T> G, A, B, C, D, E, HI, Gq, GBj, RHO;
  mat_ZZ Bmat, junk;
  vec_ZZ Rvec, Qvec;

  vec_long Rideals;
  long numRideals = 0;
  QuadraticInfElement<T, U> F, RHOdist;
  U sqReg, Fstep, curr_dist, diff;
  bool done;

  HashEntryInt<T, U> *Inode;
  IndexedHashTable<HashEntryInt<T, U>> QT, RT;
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

    if (-0.000000001 < regulator - (sqReg + get_dist_mod(delta)))
      nuclose(F, sqReg);

    if (F.is_one())
      Fstep = regulator;
    else
      Fstep = F.get_distance();

    if (DBG_CGBGRL)
      std::cout << "sqReg = " << sqReg << std::endl;
    if (DBG_CGBGRL)
      std::cout << "regulator = " << regulator << std::endl;
    if (DBG_CGBGRL)
      std::cout << "mod = " << get_dist_mod(delta) << std::endl;
    if (DBG_CGBGRL)
      std::cout << "Fstep = " << Fstep << std::endl;
    if (DBG_CGBGRL)
      std::cout << "F = " << F << std::endl;

    // initialize sets R and Q
    temp = SqrRoot(hstar << 1);
    upper = to<long>(temp);

    QT.initialize(upper << 1);
    B.assign_one();
    QT.hash(B.hash_int(to<U>(0)));

    RT.initialize(upper << 1);
    //    RT.hash (B.hash_int (to<S>(0)));

    Rideals[0] = 0;
    numRideals = 1;

    RHO.assign(B);
    RHOdist.assign(B);
    do {
      RT.hash(RHO.hash_int(to<U>(0)));
      RHOdist.rho();
      RHO.assign(RHOdist);
      if (F.is_one())
        done = (RHO.is_one());
      else
        done = (RHOdist.get_distance() > sqReg + diff || RHO.is_one());
    } while (!done);
  }

  fact_base = new QuadraticClassGroupElement<T>[100];
  contributors = new long[100];

  num_prime_ideals = 0;
  crank = 0;
  ::set(det);
  ::set(hS);
  G.assign_one();

  // while h < hstar, we only have a subgroup
  while (hS < hstar) {
    // get prime ideal
    do {
      get_next_prime(G);
      ++num_prime_ideals;
    } while (is_principal(G));

    fact_base[crank] = G;
    contributors[crank] = crank;

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

    // check whether current ideal is in previously generated subgroup
    for (i = 0; i < numQ && ::IsZero(Bjj); ++i) {
      D.assign(QT[i]);
      nucomp_real(E, D, G);

      // test for all RHO equivalent to E
      if (F.is_one()) {
        // entire equivalence classes are stored in RT
        Inode = RT.search(E.hash_int(to<U>(0)));
        if (Inode) {
          ::set(Bjj);
        }
      } else {
        // use baby-step giant-step for equivalence testing
        // baby steps are all in RT
        RHOdist.assign(E);
        ::clear(curr_dist);
        do {
          RHO.assign(RHOdist);
          Inode = RT.search(RHO.hash_int(to<U>(0)));
          if (Inode) {
            ::set(Bjj);
          }

          if (::IsZero(Bjj)) {
            curr_dist += Fstep;
            nucomp(RHOdist, RHOdist, F);
            RHOdist.adjust(curr_dist);
          }
        } while (RHOdist.get_distance().eval() <= regulator && ::IsZero(Bjj));
      }
    }

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
        nucomp_real(A, A, HI);

        if (DBG_CGBGRL)
          std::cout << "\nr = " << r << ", A = " << A << std::endl;

        for (k = 0; k < numRideals; ++k) {
          D.assign(RT[Rideals[k]]);
          nucomp_real(E, D, A);

          if (curr_index == Rideals.MaxLength())
            Rideals.SetLength(Rideals.MaxLength() << 1);

          Rideals[to<long>(curr_index)] = RT.no_of_elements();

          if (DBG_CGBGRL)
            std::cout << "BABY (NS):  r = " << r << ", E = " << E << std::endl;

          RHO.assign(E);
          RHOdist.assign(E);
          do {
            RT.hash(RHO.hash_int(curr_index));
            RHOdist.rho();
            RHO.assign(RHOdist);
            if (F.is_one())
              done = (RHO == E);
            else
              done = (RHOdist.get_distance().eval() > sqReg + diff || RHO == E);
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
          nucomp_real(D, E, B);

          // test for all RHO equivalent to E
          if (F.is_one()) {
            // entire equivalence classes are stored in RT
            Inode = RT.search(D.hash_int(to<U>(0)));
            if (Inode && !D.is_one()) {

              if (DBG_CGBGRL)
                std::cout << "FOUND GIANT" << std::endl;

              r = Inode->get_d();
              q = i;
              Bjj = y + r / numRideals;

              if (DBG_CGBGRL)
                std::cout << "r = " << r << ", q = " << q << ", Bjj = " << Bjj
                          << std::endl;

              if (Bjj > 1) {
                r %= to<U>(numRideals);
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
              RHO.assign(RHOdist);
              Inode = RT.search(RHO.hash_int(to<U>(0)));
              if (Inode && !D.is_one()) {

                if (DBG_CGBGRL)
                  std::cout << "FOUND GIANT" << std::endl;

                r = Inode->get_d();
                q = i;
                Bjj = y + r / numRideals;

                if (DBG_CGBGRL)
                  std::cout << "r = " << r << ", q = " << q << ", Bjj = " << Bjj
                            << std::endl;

                if (Bjj > 1) {
                  r %= to<U>(numRideals);
                  decode_vector(Bmat, to<ZZ>(Bjj), to<ZZ>(r), to<ZZ>(q), Rvec,
                                Qvec, numRideals, numQ);
                  // PART_OF_OLD_DEBUG: test_relation_cg(fact_base, crank,
                  // Bmat);
                }
              }

              if (::IsZero(Bjj)) {
                curr_dist += Fstep;
                nucomp(RHOdist, RHOdist, F);
                RHOdist.adjust(curr_dist);
              }
            } while (RHOdist.get_distance() <= regulator &&
                     ::IsZero(Bjj));
          }
        }

        if (::IsZero(Bjj)) {
          // not found, take another giant step
          y += u;
          nucomp_real(B, B, C);
        }
      }

      if (DBG_CGBGRL)
        std::cout << "Got giant" << std::endl;

      // double u
      if (::IsZero(Bjj)) {
        s = u + 1;
        u <<= 1;
        nudupl_real(C, C);
      }
    }

    if (!::IsOne(Bjj)) {

      hS *= Bjj;
      //	  nFI /= to<RR> (Bjj);
      ++crank;
    }

    if (DBG_CGBGRL)
      std::cout << "h = " << hS << std::endl;

    if (hS < hstar) {
      Bj = to<U>(CeilToZZ(sqrt(to<RR>(Bjj))));

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
            nucomp_real(D, E, Gq);
            QT.hash(D.hash_int(curr_index));
            ++curr_index;
          }
          nucomp_real(Gq, Gq, GBj);
          ++Qvec[crank - 1];
        }
        if (DBG_CGBGRL)
          std::cout << "DONE Q" << std::endl;
      }
    }
  }

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
      Bjj = to<U>(SNF[j][j]);
      if (!::IsOne(Bjj)) {
        CL.SetLength(i + 1);
        CL[i++] = to<ZZ>(Bjj);
      }
    }
  }
}

} // namespace ANTL
#endif
