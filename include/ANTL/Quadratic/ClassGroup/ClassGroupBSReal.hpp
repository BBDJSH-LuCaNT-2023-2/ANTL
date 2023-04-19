#ifndef CLASSGROUP_BS_REAL_H
#define CLASSGROUP_BS_REAL_H

#include <ANTL/common.hpp>
#include <unordered_set>

#include <ANTL/HashTable/HashEntryInt.hpp>
#include <ANTL/HashTable/IndexedHashTable.hpp>

#include <ANTL/LinearAlgebra/Smith.hpp>

#include <ANTL/Quadratic/QuadraticClassGroupElement.hpp>
#include <ANTL/Quadratic/QuadraticInfElement.hpp>

#include "SequencedUnorderedMap/SequencedUnorderedMap.hpp"

// template <> struct std::hash<QuadraticIdealBase<ZZ>> {
//   std::size_t operator()(QuadraticIdealBase<ZZ> const &qib) const noexcept {
//     std::size_t h1 = std::hash<int>{}(to<int>(qib.get_a()));
//     std::size_t h2 = std::hash<int>{}(to<int>(qib.get_b()));
//     return h1 ^ (h2 << 1); // or use boost::hash_combine
//   }
// };

NTL_CLIENT;
using namespace ANTL;

namespace ANTL {

// Partial class specializtion as a temporary work around multi-pararameter
// template restrictions.
template <class T> class ClassGroupBSReal {
private:
  // DBG_CONSTANTS
  bool DBG_CGBSRL = false;
  bool DBG_CGBSRL_SUM = false;

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
  vector<ZZ> CL;

  // size of fact_base
  long numFB;

  // transformation matrix
  mat_ZZ U;

  bool reset_prime_seq;

public:
  ClassGroupBSReal(QuadraticOrder<T> *quadratic_order);

  ~ClassGroupBSReal();

  void cg_bs_real(const ZZ &hstar);

  void set_regulator(double &ext_regulator){regulator = ext_regulator;}

  vector<ZZ> get_class_group(){return CL;}

private:
  ZZ get_dist_mod(const T &Delta) { return CeilToZZ(log(to_RR(to_ZZ(Delta)))); }

  void get_next_prime(QuadraticClassGroupElement<T> &G);
};

// Method definitions - Everything below will eventually go into a
// RegulatorLenstraData_impl.hpp file.

template <class T>
ClassGroupBSReal<T>::ClassGroupBSReal(QuadraticOrder<T> *quadratic_order_arg) {

  quadratic_order = quadratic_order_arg;
  delta = quadratic_order->get_discriminant();

  //   parallel = false;
  //   Rbsgs = false;        // true if R was computed using BSGS
  //   Rconditional = false; // true if correctness of R relies on ERH
}

template <class T>
ClassGroupBSReal<T>::~ClassGroupBSReal() {
  delete[] fact_base;
  delete[] contributors;
}

template <class T> void ClassGroupBSReal<T>::cg_bs_real(const ZZ &hstar) {
  if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 0 - BEGIN" << std::endl;}
  long i, j, k, l, n, m, x, q, r, upper, size1, size2, index, hashNum, curr;
  long sizeH1, sizeH2, sizeAuxBabySet, H1Index, H2Index, auxIndex;
  long rem, divide, size;
  ZZ s, t, e, temp, product1, product2, vIndex, wIndex, Bjj;
  bool found;

  reset_prime_seq = true;

  QuadraticClassGroupElement<T> giantElement{*quadratic_order}, babyElement{*quadratic_order}, RHO{*quadratic_order};
  QuadraticClassGroupElement<T> g{*quadratic_order}, gjx{*quadratic_order}, ginv{*quadratic_order}, h1{*quadratic_order}, h2{*quadratic_order}, a{*quadratic_order}, a_test{*quadratic_order}, b{*quadratic_order}, c{*quadratic_order};

  mat_ZZ B, tempB, SNF;
  vec_ZZ v, w, Bj, H1Vec, H2Vec;

  vec_long Rideals;
  long numRideals;
  QuadraticInfElement<T, double> F{*quadratic_order}, RHOdist{*quadratic_order};
  ZZ sqReg, Fstep, curr_dist, diff;
  bool done;

  sizeH1 = sizeH2 = sizeAuxBabySet = numRideals = 0;

  HashEntryInt<T, ZZ> *node, *node_test;
//   IndexedHashTable<HashEntryInt<T, ZZ>> ibabySet, igiantSet;
  SequencedUnorderedMap<T> ibabySet, igiantSet;

  long *I1, *I2;
  if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 1 - Initialize variables" << std::endl;}
  // initialize variables
  j = 0;
  l = 100;
  m = -1;
  size1 = 0;
  size2 = 0;
  num_prime_ideals = 0;
  set(h);
  clear(s);
  clear(t);
  found = false;

  if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 2 - initialize data structures" << std::endl;}
  // initialize data structures
  B.SetDims(0, 0);
  fact_base = new QuadraticClassGroupElement<T>[l];
  contributors = new long[l];
  I1 = new long[l];
  I2 = new long[l];

  if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3 - initialize hash tables H1 and H2" << std::endl;}
  // initialize hash tables H1 and H2
  temp = SqrRoot(hstar << 1);
  conv(upper, temp);

  if (hstar > 1) {
    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.0" << std::endl;}
    // initialize prime ideal
    g.assign_one();

    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.1" << std::endl;}
    // initialize hash tables H1 and H2
    temp = SqrRoot(hstar << 1);
    conv(upper, temp);

    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.2" << std::endl;}
//     ibabySet.initialize(upper << 1);
//     igiantSet.initialize(upper << 1);

    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.3" << std::endl;}
    a.assign_one();
    sizeH1 = sizeH2 = sizeAuxBabySet = 1;
    igiantSet.hash(a.hash_int(ZZ::zero()));

    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.4.0" << std::endl;}
    // get giant step for equivalence testing
    SqrRoot(sqReg, regulator);
    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.4.1" << std::endl;}
    diff = get_dist_mod(delta);
    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.4.2" << std::endl;}
    F.assign_one();

    if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 3.5" << std::endl;}
    if (-0.000000001 < regulator - to<double>(sqReg + diff)) {
      nuclose(F, sqReg);
    }

    if (F.is_one()) {
      Fstep = regulator;
    }
    else {
      Fstep = F.get_distance();
    }

    Rideals.SetLength(1024);
    Rideals[0] = 0;
    numRideals = 1;

    RHO.assign(a);
    RHOdist.assign(a);

    do {
      ibabySet.hash(RHO.hash_int(ZZ::zero()));
//       ibabySet_test.hash(RHO.hash_int(ZZ::zero()));
//       RHOdist.rho();
      RHOdist.baby_step();
      RHO.assign(RHOdist.get_qib());

      if (F.is_one())
        done = (RHO.IsOne());
      else
//         done = (RHOdist.get_distance() > (sqReg + diff) || RHO.IsOne());
        done = (0.000000001 > to<double>(sqReg + diff) - RHOdist.get_distance() || RHO.IsOne());
    } while (!done);
  }

  if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 4 - iterate until h > hstar, otherwise we only have a subgroup" << std::endl;}
  // iterate until h > hstar, otherwise we only have a subgroup
  while (h < hstar) {
#ifdef DEBUGBS
    cout << "p" << flush;
#endif
    // get a prime ideal
    get_next_prime(g);
    num_prime_ideals++;
    fact_base[j] = g;
    contributors[j] = j;

    // initialize elements
    set(e);
    babyElement = g;
    giantElement = g;
    found = false;

    // loop until we find a match
    while (!found) {
      // for all (g,v) in giantSet
      hashNum = igiantSet.no_of_elements();
      if (DBG_CGBSRL_SUM) {
            std::cout << "hashNum is " << hashNum << std::endl;
      }

#ifdef DEBUGBS
      cout << "2a" << flush;
#endif
      for (i = 0; i < hashNum; i++) {
#ifdef DEBUGBS
        cout << "2b" << flush;
#endif
        a.assign(igiantSet[i]);
        if (DBG_CGBSRL_SUM) {
          std::cout << "a is " << a << std::endl;
        }
//         nucomp_real(b, giantElement, a);
        mul(b, giantElement, a);
        vIndex = igiantSet[i].get_d();
        if (DBG_CGBSRL_SUM) {
          std::cout << "vIndex is " << vIndex << std::endl;
        }

        // test for all RHO equivalent to b
        if (F.is_one()) {
          if (DBG_CGBSRL_SUM) {
            ibabySet.print(std::cout);
            std::cout << "looking for" << b.hash_int(ZZ::zero());
          }
          node = ibabySet.search(b.hash_int(ZZ::zero()));
          if (DBG_CGBSRL_SUM) {
            std::cout << "node is " << node << std::endl;
          }
//           node_test = ibabySet_test.search(b.hash_int(ZZ::zero()));
        }
        else {
          // use baby-step giant-step for equivalence testing
          // baby steps are all in ibabySet
          RHOdist.assign(b);
          clear(curr_dist);
          do {
            /*		  #ifdef DEBUGBS
              cout << RHOdist.get_distance().eval() << " up to " << R.eval() <<
              endl; #endif*/
            RHO.assign(RHOdist.get_qib());
            node = ibabySet.search(RHO.hash_int(ZZ::zero()));
            if (DBG_CGBSRL_SUM) {
              std::cout << "node is " << node << std::endl;
            }
//             node_test = ibabySet_test.search(RHO.hash_int(ZZ::zero()));

            if (node)
              break;

            curr_dist += Fstep;
//             nucomp(RHOdist, RHOdist, F);
            RHOdist.giant_step(F);
            RHOdist.adjust(curr_dist);
          } while (-0.000000001 <= regulator - RHOdist.get_distance());
        }

        if (node) {
          // get w vector
          wIndex = node->get_d();
          w.SetLength(j + 1);
          clear(w);

          auxIndex = wIndex % sizeAuxBabySet;
          H1Index = auxIndex % sizeH1;

          Bjj = e * (e + 1) / 2 - (wIndex / sizeAuxBabySet);

          if (Bjj > 1) {
            // find corresponding H1 vector
            rem = H1Index;
            size = H1Vec.length();
            divide = sizeH1;
            for (k = j - 1; k >= 0; k--) {
              if (rem == 0) {
                clear(w[k]);
              } else {
                divide /= to_int(H1Vec[k]);
                w[k] = rem / divide;
                rem %= divide;
              }
            }

            // find corresponding auxBabySet vector
            if (m >= 0) {
              w[m] += (auxIndex / sizeH1);
            }

            // find corresponding babySet vector
            w[j] -= (wIndex / sizeAuxBabySet);

            // get v vector
            v.SetLength(j + 1);
            clear(v);

            H2Index = vIndex % sizeH2;

            // find corresponding H2 vector
            rem = H2Index;
            divide = sizeH2;
            for (k = j - 1; k >= 0; k--) {
              if (rem == 0) {
                clear(v[k]);
              } else {
                divide /= to_int(H2Vec[k]);
                v[k] = rem / divide;
                rem %= divide;
              }
            }

            // find corresponding giantSet vector
            if (m >= 0) {
              v[m] += ((vIndex / sizeH2) * s);
            }

            // compute Bj
            Bj.SetLength(j + 1);
            Bj = v + w;
            Bj[j] += e * (e + 1) / 2;

            // resize basis matrix
            tempB.SetDims(j + 1, j + 1);
            for (k = 0; k < j; k++) {
              for (n = 0; n < j; n++) {
                tempB[k][n] = B[k][n];
              }
            }

            // copy over new vector
            for (k = 0; k <= j; k++) {
              tempB[k][j] = Bj[k];
            }

            // reset basis matrix
            B.kill();
            B.SetDims(j + 1, j + 1);
            B = tempB;
          } // end if (Bjj > 1)

          found = true;
          break;

        } // end if (node)
      }

      if (found)
        break;
#ifdef DEBUGBS
      cout << "CANT_FIND" << flush;
#endif

      // update ibabySet
      curr = numRideals;

      for (i = 0; i < sizeAuxBabySet; i++) {
#ifdef DEBUGBS
        cout << "3" << flush;
#endif
        a.assign(ibabySet[Rideals[i]]);
        if (DBG_CGBSRL_SUM) {
          std::cout << "a is " << a << std::endl;
        }
//         a_test.assign(ibabySet_test.at(Rideals[i]));
//         nucomp_real(b, a, babyElement);
        mul(b, a, babyElement);

        if (curr == Rideals.MaxLength())
          Rideals.SetLength(Rideals.MaxLength() << 1);
        if(DBG_CGBSRL_SUM) {
          std::cout << "ibabySet.no_of_elements() is " << ibabySet.no_of_elements() <<std::endl;
        }
        Rideals[to<long>(curr)] = ibabySet.no_of_elements();

        RHO.assign(b);
        RHOdist.assign(b);
        do {
          ibabySet.hash(RHO.hash_int(to_ZZ(curr)));
//           ibabySet_test.hash(RHO.hash_int(to_ZZ(curr)));
          RHOdist.baby_step();
          RHO.assign(RHOdist.get_qib());

          if (F.is_one())
            done = (RHO == b);
          else
            done = (RHOdist.get_distance() > sqReg + diff || RHO == b);
        } while (!done);
        curr++;
      }

      // update variables
      e++;
//       nucomp_real(babyElement, babyElement, g);
      mul(babyElement, babyElement, g);
//       nucomp_real(giantElement, giantElement, babyElement);
      mul(giantElement, giantElement, babyElement);
      numRideals = curr;

    } // end while (!found)

#ifdef DEBUGBS
    cout << "d" << flush;
#endif
    // update h
    h *= Bjj; //[j][j];

    if (h < hstar && Bjj > 1) {
      // delete all elements of ibabySet until only H1 remains
      if (sizeH1 < numRideals) {
        hashNum = ibabySet.no_of_elements() - 1;
        if(DBG_CGBSRL_SUM) {
          std::cout << "hashNum is " << hashNum <<std::endl;
        }
        for (i = hashNum; i >= Rideals[sizeH1]; i--) {
          ibabySet.remove_from(i);
//           ibabySet_test.remove_from(i);
        }
      }

      // delete all elements of igiantSet until only H2 remains
      hashNum = igiantSet.no_of_elements();
      if(DBG_CGBSRL_SUM) {
        std::cout << "hashNum is " << hashNum << std::endl;
      }
      if (sizeH2 < hashNum) {
        for (i = hashNum - 1; i >= sizeH2; i--) {
          igiantSet.remove_from(i);
        }
      }

      // update H1 and H2 vectors
      H1Vec.SetLength(j + 1);
      H2Vec.SetLength(j + 1);
      H1Vec[j] = 1;
      H2Vec[j] = 1;

      // compute product(Bii) for I1 and I2
      product1 = 1;
      for (i = 0; i < size1; i++) {
        index = I1[i];
        product1 *= B[index][index];
      }

      product2 = 1;
      for (i = 0; i < size2; i++) {
        index = I2[i];
        product2 *= B[index][index];
      }

      // decide which set to add elements to
      if (B[j][j] * product1 <= SqrRoot(h)) {
        // add to set H1
        curr = sizeH1;
        ginv.assign(conjugate(g));
        gjx.assign_one();

        for (x = 1; x < B[j][j]; x++) {
//           nucomp_real(gjx, gjx, ginv);
          mul(gjx, gjx, ginv);
#ifdef DEBUGBS
          cout << "5a0" << flush;
#endif

          for (i = 0; i < sizeH1; i++) {
#ifdef DEBUGBS
            cout << "5a1" << flush;
#endif
            h1.assign(ibabySet[Rideals[i]]);
//             nucomp_real(b, gjx, h1);
            mul(b, gjx, h1);

            if (curr == Rideals.MaxLength())
              Rideals.SetLength(Rideals.MaxLength() << 1);
            if(DBG_CGBSRL_SUM) {
              std::cout << "ibabySet.no_of_elements() is " << ibabySet.no_of_elements() << std::endl;
            }
            Rideals[to<long>(curr)] = ibabySet.no_of_elements();

            RHO.assign(b);
            RHOdist.assign(b);
            do {
              ibabySet.hash(RHO.hash_int(to_ZZ(curr)));
              RHOdist.baby_step();
              RHO.assign(RHOdist.get_qib());

              if (F.is_one())
                done = (RHO == b);
              else
                done =
                    (RHOdist.get_distance() > sqReg + diff || RHO == b);
            } while (!done);
            curr++;
          }
        }
        sizeH1 = curr;

        // update indices
        I1[size1] = j;
        size1++;
        product1 *= B[j][j];
        H1Vec[j] = B[j][j];

      } else {
        // add to set H2
        if (m >= 0) {
          curr = sizeH2;
          gjx.assign_one();

          for (x = 1; x < B[m][m]; x++) {
//             nucomp_real(gjx, gjx, fact_base[m]);
            mul(gjx, gjx, fact_base[m]);
#ifdef DEBUGBS
            cout << "5b0" << flush;
#endif
            for (i = 0; i < sizeH2; i++) {
#ifdef DEBUGBS
              cout << "5b1" << flush;
#endif
              h2.assign(igiantSet[i]);
//               nucomp_real(b, gjx, h2);
              mul(b, gjx, h2);
              igiantSet.hash(b.hash_int(to_ZZ(curr)));
              curr++;
            }
          }
          sizeH2 = curr;

          // update indices
          I2[size2] = m;
          size2++;
          product2 *= B[m][m];
          H2Vec[m] = B[m][m];
          H1Vec[m] = 1;
        }

        m = j;
      }

      // update s and t
      s = CeilToZZ(sqrt(to_RR(h)) / to_RR(product1));
      t = CeilToZZ(sqrt(to_RR(h)) / to_RR(product2));

      // update ibabySet
      curr = sizeH1;
      ginv.assign(conjugate(fact_base[m]));
      gjx.assign_one();

      for (r = 1; r < s; r++) {
#ifdef DEBUGBS
        cout << "7a0" << flush;
#endif
//         nucomp_real(gjx, gjx, ginv);
        mul(gjx, gjx, ginv);
        for (i = 0; i < sizeH1; i++) {
#ifdef DEBUGBS
          cout << "7a1" << flush;
#endif
          h1.assign(ibabySet[Rideals[i]]);
//           nucomp_real(b, h1, gjx);
          mul(b, h1, gjx);

          if (curr == Rideals.MaxLength())
            Rideals.SetLength(Rideals.MaxLength() << 1);
          if(DBG_CGBSRL_SUM) {
              std::cout << "ibabySet.no_of_elements() is " << ibabySet.no_of_elements() << std::endl;
          }
          Rideals[to<long>(curr)] = ibabySet.no_of_elements();

          RHO.assign(b);
          RHOdist.assign(b);
          do {
            ibabySet.hash(RHO.hash_int(to_ZZ(curr)));
            RHOdist.baby_step();
            RHO.assign(RHOdist.get_qib());

            if (F.is_one())
              done = (RHO == b);
            else
              done = (0.000000001 > to<double>(sqReg + diff) - RHOdist.get_distance() || RHO == b);

          } while (!done);
          curr++;
        }
      }

      sizeAuxBabySet = numRideals = curr;

      // update igiantSet
      curr = sizeH2;
      nupower_real(c, fact_base[m], s);
      gjx.assign_one();

      for (q = 1; q < t; q++) {
//         nucomp_real(gjx, gjx, c);
        mul(gjx, gjx, c);
        for (i = 0; i < sizeH2; i++) {
          h2.assign(igiantSet[i]);
          if(DBG_CGBSRL_SUM) {
            std::cout << "h2 is " << h2 << std::endl;
          }
//           nucomp_real(b, h2, gjx);
          mul(b, h2, gjx);
          igiantSet.hash(b.hash_int(to_ZZ(curr)));
          curr++;
        }
      }
    } // end if (Bjj > 1)

    // if B[j][j] > 1 we store the generator and move on
    if (Bjj > 1) {
      j++;
    }

    // check if data structures are full
    if (j == l) {
      QuadraticClassGroupElement<T> *temp1;
      long *temp2, *tempI1, *tempI2;

      temp1 = fact_base;
      temp2 = contributors;
      tempI1 = I1;
      tempI2 = I2;
      fact_base = new QuadraticClassGroupElement<T>[j << 1];
      contributors = new long[j << 1];
      I1 = new long[j << 1];
      I2 = new long[j << 1];

      for (i = 0; i < j; i++) {
        fact_base[i] = temp1[i];
        contributors[i] = temp2[i];
        I1[i] = tempI1[i];
        I2[i] = tempI2[i];
      }

      delete[] temp1;
      delete[] temp2;
      delete[] tempI1;
      delete[] tempI2;

      l = j << 1;
    }

  } // end while (h <= hstar)

  if (DBG_CGBSRL) {std::cout << "CGBSRL: STEP 5 - compute the structure" << std::endl;}
  // compute the structure
  if (h == 1) {
    CL.push_back(ZZ(1));
//     fact_base[0].assign_one();
    QuadraticClassGroupElement<T> identity{*quadratic_order};
    fact_base[0] = identity;
    contributors[0] = 0;
    numFB = 1;
  } else {
    numFB = j;
    SNF = SmithNormalForm(B, U, tempB);

    i = 0;
    for (j = 0; j < numFB; j++) {
      temp = SNF[j][j];
      if (!IsOne(temp)) {
        CL.push_back(ZZ(temp));
      }
    }
  }

  delete[] I1;
  delete[] I2;
}

template <class T>
void ClassGroupBSReal<T>::get_next_prime(QuadraticClassGroupElement<T> &G) {
  static PrimeSeq PS;

//   if (G.IsOne()) {
//     PS.reset(0);
//   }

    if (reset_prime_seq) {
//     std::cout << "resetting PS!" << std::endl;
    PS.reset(0);
    reset_prime_seq = false;
  }

  long p = PS.next();

  while (!G.assign_prime(to<T>(p)))
    p = PS.next();
}

} // namespace ANTL
#endif
