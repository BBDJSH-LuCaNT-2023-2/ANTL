#ifndef AUXILLARYFUNCTIONS_HPP_INCLUDED
#define AUXILLARYFUNCTIONS_HPP_INCLUDED

#define INT_SIZE_QUADRATIC 10000000

#include <algorithm>
#include <iostream>
#include <list>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#include <ANTL/Quadratic/ClassGroup/ClassGroupBSGSReal.hpp>
#include "../../../../src/Quadratic/ClassGroup/ClassGroupBSGSReal.cpp"

#include <ANTL/Quadratic/ClassGroup/ClassGroupBSReal.hpp>

#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_long.hpp>

#include <ANTL/Quadratic/QuadraticClassGroupElement.hpp>

using namespace NTL;
using namespace ANTL;

pair<double, ZZ> get_regulator_and_hstar(QuadraticOrder<long> &quad_order, L_function<long> &l_function) {
  // Setting up the RegulatorLenstraData object
  RegulatorLenstraData<long, double> regulator_lenstra_data{&quad_order,
                                                          &l_function};

  // Computing the regulator
  regulator_lenstra_data.regulator_lenstra();

  double regulator = regulator_lenstra_data.get_regulator();

  // Computing h*
  RR h_star_close = to_RR(regulator_lenstra_data.lower_bound_hR()) / to_RR(regulator);
  ZZ h_star = CeilToZZ(h_star_close);

  return {regulator, h_star};
}

vector<long> get_class_group_BSGS(QuadraticOrder<long> &quad_order, double regulator, ZZ h_star) {
  // Setting up the ClassGroupBSGSReal object
  ClassGroupBSGSReal<long> class_group_bsgs_real1{&quad_order};
  class_group_bsgs_real1.set_regulator(regulator);

  // Computing the class group
  class_group_bsgs_real1.cg_bsgs_real(h_star);

  // Adding computed class group to reslults vector
  vector<ZZ> class_group_ZZ = class_group_bsgs_real1.get_class_group();
  vector<long> class_group_long = {};

  for(auto num : class_group_ZZ) {
    class_group_long.push_back(to<long>(num));
  }

  std::sort(class_group_long.begin(), class_group_long.end());

  return class_group_long;
}

vector<long> get_class_group_BS(QuadraticOrder<long> &quad_order, double regulator, ZZ h_star) {
  // Setting up the ClassGroupBSGSReal object
  ClassGroupBSReal<long> class_group_bs_real1{&quad_order};
  class_group_bs_real1.set_regulator(regulator);

  // Computing the class group
  class_group_bs_real1.cg_bs_real(h_star);

  // Adding computed class group to reslults vector
  vector<ZZ> class_group_ZZ = class_group_bs_real1.get_class_group();
  vector<long> class_group_long = {};

  for(auto num : class_group_ZZ) {
    class_group_long.push_back(to<long>(num));
  }

  std::sort(class_group_long.begin(), class_group_long.end());

  return class_group_long;
}

void get_DList_real_custom(long ubound, std::list<long> &discriminants);

void generateDiscriminantsConcurrent(
    long ubound, int offset, int num_of_threads,
    std::unordered_set<long> &local_discriminants);
bool isDiscriminant(long &candidate, std::list<long> &primes);
bool isSquareFree(long &candidate, std::list<long> &primes);
void generatePrimes(std::list<long> &primes, long &bound);

void generateDiscriminantsConcurrent(
    long ubound, int offset, int thread_offset,
    std::unordered_set<long> &local_discriminants) {

  // Generate fundamental discriminants
  std::list<long> local_primes;
  long bound = ceil(sqrt(ubound));
  generatePrimes(local_primes, bound);

  std::vector<long> thread_candidates(6);
  thread_candidates.at(0) = offset + 1;
  thread_candidates.at(1) = offset + 5;
  thread_candidates.at(2) = offset + 8;
  thread_candidates.at(3) = offset + 9;
  thread_candidates.at(4) = offset + 12;
  thread_candidates.at(5) = offset + 13;

  while (thread_candidates.at(5) < ubound) {
    if (isDiscriminant(thread_candidates.at(0), local_primes)) {
      local_discriminants.insert(thread_candidates.at(0));
    }
    if (isDiscriminant(thread_candidates.at(1), local_primes)) {
      local_discriminants.insert(thread_candidates.at(1));
    }
    if (isDiscriminant(thread_candidates.at(2), local_primes)) {
      local_discriminants.insert(thread_candidates.at(2));
    }
    if (isDiscriminant(thread_candidates.at(3), local_primes)) {
      local_discriminants.insert(thread_candidates.at(3));
    }
    if (isDiscriminant(thread_candidates.at(4), local_primes)) {
      local_discriminants.insert(thread_candidates.at(4));
    }
    if (isDiscriminant(thread_candidates.at(5), local_primes)) {
      local_discriminants.insert(thread_candidates.at(5));
    }

    for (int i = 0; i < 6; i++) {
      thread_candidates.at(i) += thread_offset;
    }
  }

  for (auto thread_candidate : thread_candidates) {
    if (thread_candidate < ubound) {
      if (isDiscriminant(thread_candidate, local_primes)) {
        local_discriminants.insert(thread_candidate);
      }
    }
  }
}

bool isDiscriminant(long &candidate, std::list<long> &primes) {
  // A fundamental discriminant mod 16 == 1, 5, 8, 9, 12, 13
  long long candidateMod16 = candidate % 16;
  if (!(candidateMod16 == 1 || candidateMod16 == 5 || candidateMod16 == 8 ||
        candidateMod16 == 9 || candidateMod16 == 12 || candidateMod16 == 13)) {
    return false;
  }

  // For a fundamental discriminant, delta, it must be the case
  // thar either delta or delta/4 is square-free
  long squareCheckCandidate = candidate;

  if (candidate % 4 == 0) {
    squareCheckCandidate /= 4;
  }

  if (!isSquareFree(squareCheckCandidate, primes)) {
    return false;
  }

  // Both criteria are met
  return true;
}

bool isSquareFree(long &candidate, std::list<long> &primes) {
  // Square-free test
  for (std::list<long>::iterator i = primes.begin(); i != primes.end(); i++) {
    long long squared_primed = (*i) * (*i);
    if (squared_primed > candidate) {
      break;
    }
    if (candidate % squared_primed == 0) {
      return false;
    }
  }
  return true;
}

void generatePrimes(std::list<long> &primes, long &bound) {

  std::vector<bool> is_a_prime(bound + 1, true);

  is_a_prime[0] = false;
  is_a_prime[1] = false;

  for (long current_prime = 2; current_prime * current_prime <= bound;
       current_prime++) {
    if (is_a_prime[current_prime]) {
      for (long not_a_prime = current_prime * current_prime;
           not_a_prime <= bound; not_a_prime += current_prime) {
        is_a_prime.at(not_a_prime) = 0;
      }
    }
  }

  for (long i = 0; i < bound + 1; i++) {
    if (is_a_prime[i]) {
      primes.push_back(i);
    }
  }
}

void get_DList_real_custom(long ubound, std::list<long> &discriminants) {

  long bound = ceil(sqrt(ubound));
  std::vector<bool> is_a_prime(bound + 1, true);
  std::list<long> primes;

  is_a_prime[0] = false;
  is_a_prime[1] = false;

  for (long current_prime = 2; current_prime * current_prime <= bound;
       current_prime++) {
    if (is_a_prime[current_prime]) {
      for (long not_a_prime = current_prime * current_prime;
           not_a_prime <= bound; not_a_prime += current_prime) {
        is_a_prime.at(not_a_prime) = 0;
      }
    }
  }

  for (long i = 0; i < bound + 1; i++) {
    if (is_a_prime[i]) {
      primes.push_back(i);
    }
  }
  std::vector<long> offsets = {2, 3, 4, 6, 7, 10, 11, 14, 15};
  std::vector<bool> is_a_discriminant(ubound+1, true);
  is_a_discriminant[0] = 0;
  is_a_discriminant[1] = 0;

  for(int i = 0; i < is_a_discriminant.size(); i += 16) {
    is_a_discriminant[i] = 0;
    for(auto offset : offsets) {
      if(i+offset < is_a_discriminant.size()) {
        is_a_discriminant[i+offset] = 0;
      }
    }
  }

  for(auto prime : primes) {
    std::cout << "current prime is " << prime << std::endl;
    if(prime != 2) {
      for(long prime_square = prime * prime; prime_square < is_a_discriminant.size(); prime_square += prime * prime){
        std::cout << "current prime square is " << prime_square << std::endl;
        is_a_discriminant[prime_square] = 0;
      }
    }
  }


  for (long i = 0; i < ubound + 1; i++) {
    if (is_a_discriminant[i] == 1) {
      discriminants.push_back(i);
    }
  }

}

//
// Compute a list of real quadratic discriminants with L <= D <= H
// Assumes INT_SIZE is <= 10^9
//

void get_Dlist_real(const ZZ &L, const ZZ &H, std::vector<ZZ> &Dlist, long &n, int init) {
  static long prime_list[1049000]; // first 1000000 primes
  static long numP = 0;
  std::vector<unsigned char> sieve0(1000000000), sieve1(1000000000);
  unsigned char *end0, *end1;
  ZZ B, PP0, PP1, P2, LL0, LL1, off;
  long p, p2l, off0, off1;
  unsigned char *sptr0, *sptr1;
  long *pl;

  std::cout << "hello 1" << std::endl;
  if (init) {
    // initialize static list of small odd primes <= sqrt(H)
    PrimeSeq sP;
    long bound = to<long>(SqrRoot(H)) + 100;
    numP = 0;

    sP.reset(3);
    p = sP.next();
    while (p <= bound) {
      prime_list[numP] = p;
      ++numP;
      p = sP.next();
      if (numP > 1049000) {
        cerr << "ERROR in get_Dlist_real:  prime list too small, numP = "
             << numP << endl;
        exit(1);
      }
    }
    prime_list[numP] = p;
    ++numP;

    return;
  }

  if (numP == 0) {
    cerr << "ERROR in get_DList_real:  must call with init>0 first" << endl;
    exit(1);
  }

  std::cout << "hello 2" << std::endl;
  // initialize sieve
  LL0 = L;
  if (IsZero(LL0))
    ++LL0;
  while (rem(LL0, 4) != 0)
    ++LL0;

  LL1 = L;
  while (rem(LL1, 4) != 1)
    ++LL1;

  std::cout << "hello 2.1" << std::endl;
  conv(off0, 1 + ((H - LL0) >> 2));
  end0 = &(*sieve0.begin()) + off0;
  conv(off1, 1 + ((H - LL1) >> 2));
  end1 = &(*sieve1.begin()) + off1;


  std::cout << "hello 2.2" << std::endl;
//   memset(sieve0, '\0', (off0) * sizeof(unsigned char));
  for(long i = 0; i < sieve0.size(); i++) {
    sieve0.at(i) = '\0' ;
  }
  std::cout << "hello 2.3" << std::endl;

//   memset(sieve1, '\0', (off1) * sizeof(unsigned char));
  std::cout << "hello 2.4" << std::endl;
  for(long i = 0; i < sieve0.size(); i++) {
    sieve1.at(i) = '\0' ;
  }

  B = SqrRoot(H);

  if (B > prime_list[numP - 1]) {
    cerr << "ERROR get_Dlist_real:  not enough primes for H = " << H << endl;
    exit(1);
  }

  pl = prime_list;
  p = *pl;
  P2 = to<ZZ>(p) * to<ZZ>(p);

  std::cout << "hello 3" << std::endl;
  // sieve out all odd square factors
  while (p <= B) {
    if (P2.SinglePrecision()) {
      conv(p2l, P2);

      // process sieve0 (0 mod 4 integers)
      off0 = rem(LL0, p2l);
      if (off0)
        off0 = rem(to<ZZ>(InvMod(4, p2l)) * to<ZZ>(-off0), p2l);
      if (off0 < 0)
        off0 += p2l;

      sptr0 = &(*sieve0.begin()) + off0;
      while (sptr0 <= end0) {
        (*sptr0) = 1;
        sptr0 += p2l;
      }

      // process sieve1 (1 mod 4 integers)
      off1 = rem(LL1, p2l);
      if (off1)
        off1 = rem(to<ZZ>(InvMod(4, p2l)) * to<ZZ>(-off1), p2l);
      if (off1 < 0)
        off1 += p2l;

      sptr1 = &(*sieve1.begin()) + off1;
      while (sptr1 <= end1) {
        (*sptr1) = 1;
        sptr1 += p2l;
      }
    } else {
      rem(off, LL0, P2);
      if (!IsZero(off))
        rem(off, InvMod(to<ZZ>(4), P2) * to<ZZ>(-off), P2);
      if (off < 0)
        off += P2;

      if (off.SinglePrecision()) {
        conv(off0, off);
        sptr0 = &(*sieve0.begin()) + off0;
        if (sptr0 >= &(*sieve0.begin()) && sptr0 <= end0)
          (*sptr0) = 1;
      }

      // process sieve1 (1 mod 4 integers)
      rem(off, LL1, P2);
      if (!IsZero(off))
        rem(off, InvMod(to<ZZ>(4), P2) * to<ZZ>(-off), P2);
      if (off < 0)
        off += P2;

      if (off.SinglePrecision()) {
        conv(off1, off);
        sptr1 = &(*sieve1.begin()) + off1;
        if (sptr1 >= &(*sieve1.begin()) && sptr1 <= end1)
          (*sptr1) = 1;
      }
    }
    ++pl;
    p = *pl;
    P2 = to<ZZ>(p) * to<ZZ>(p);
  }

  std::cout << "hello 4" << std::endl;
  // sieve1 represents squarefree ints that are 1 mod 4
  // sieve0 represetns ints 0 mod 4 with no odd square factors.  Record those
  // for which
  //   x / 4 is 2 or 3 mod 4

  // get first D = 0 mod 4
  PP0 = LL0;
  sptr0 = &(*sieve0.begin());
  while (*sptr0) {
    PP0 += 4;
    ++sptr0;
  }

  // get first D = 1 mod 4
  PP1 = LL1;
  sptr1 = &(*sieve1.begin());
  while (*sptr1) {
    PP1 += 4;
    ++sptr1;
  }

  n = 0;
  std::cout << "hello 5" << std::endl;
  while (PP0 <= H || PP1 <= H) {
    if (PP1 < PP0 && PP1 <= H) {
      Dlist[n] = PP1;
      ++n;

      // get next D = 1 mod 4
      do {
        PP1 += 4;
        ++sptr1;
      } while (PP1 <= H && (*sptr1));
    } else {
      off = PP0 >> 2;
      p2l = rem(off, 4);
      if (p2l == 2 || p2l == 3) {
        Dlist[n] = PP0;
        ++n;
      }

      // get next D = 0 mod 4
      do {
        PP0 += 4;
        ++sptr0;
      } while (PP0 <= H && (*sptr0));
    }
  }
}

void get_Dlist_real(const ZZ &L, const ZZ &H, ZZ *Dlist, long &n, int init) {
  static long prime_list[1049000]; // first 1000000 primes
  static long numP = 0;
  static unsigned char sieve0[INT_SIZE_QUADRATIC], sieve1[INT_SIZE_QUADRATIC];
  unsigned char *end0, *end1;
  ZZ B, PP0, PP1, P2, LL0, LL1, off;
  long p, p2l, off0, off1;
  unsigned char *sptr0, *sptr1;
  long *pl;


  if (init) {
    // initialize static list of small odd primes <= sqrt(H)
    PrimeSeq sP;
    long bound = to<long>(SqrRoot(H)) + 100;
    numP = 0;

    sP.reset(3);
    p = sP.next();
    while (p <= bound) {
      prime_list[numP] = p;
      ++numP;
      p = sP.next();
      if (numP > 1049000) {
        std::cout << "get_Dlist_real init: error" << std::endl;
        cerr << "ERROR in get_Dlist_real:  prime list too small, numP = "
             << numP << endl;
        exit(1);
      }
    }
    prime_list[numP] = p;
    ++numP;
    return;
  }

  if (numP == 0) {
    cerr << "ERROR in get_DList_real:  must call with init>0 first" << endl;
    exit(1);
  }

  // initialize sieve
  LL0 = L;
  if (IsZero(LL0))
    ++LL0;
  while (rem(LL0, 4) != 0)
    ++LL0;

  LL1 = L;
  while (rem(LL1, 4) != 1)
    ++LL1;


  conv(off0, 1 + ((H - LL0) >> 2));
  end0 = sieve0 + off0;
  conv(off1, 1 + ((H - LL1) >> 2));
  end1 = sieve1 + off1;

  memset(sieve0, '\0', (off0) * sizeof(unsigned char));
  memset(sieve1, '\0', (off1) * sizeof(unsigned char));

  B = SqrRoot(H);


  if (B > prime_list[numP - 1]) {
    cerr << "ERROR get_Dlist_real:  not enough primes for H = " << H << endl;
    exit(1);
  }

  pl = prime_list;
  p = *pl;
  P2 = to<ZZ>(p) * to<ZZ>(p);

  // sieve out all odd square factors
  while (p <= B) {
    if (P2.SinglePrecision()) {
      conv(p2l, P2);



      // process sieve0 (0 mod 4 integers)
      off0 = rem(LL0, p2l);
      if (off0)
        off0 = rem(to<ZZ>(InvMod(4, p2l)) * to<ZZ>(-off0), p2l);
      if (off0 < 0)
        off0 += p2l;

      sptr0 = sieve0 + off0;
      while (sptr0 <= end0) {
        (*sptr0) = 1;
        sptr0 += p2l;
      }

      // process sieve1 (1 mod 4 integers)
      off1 = rem(LL1, p2l);
      if (off1)
        off1 = rem(to<ZZ>(InvMod(4, p2l)) * to<ZZ>(-off1), p2l);
      if (off1 < 0)
        off1 += p2l;

      sptr1 = sieve1 + off1;
      while (sptr1 <= end1) {
        (*sptr1) = 1;
        sptr1 += p2l;
      }
    } else {


      rem(off, LL0, P2);
      if (!IsZero(off))
        rem(off, InvMod(to<ZZ>(4), P2) * to<ZZ>(-off), P2);
      if (off < 0)
        off += P2;

      if (off.SinglePrecision()) {
        conv(off0, off);
        sptr0 = sieve0 + off0;
        if (sptr0 >= sieve0 && sptr0 <= end0)
          (*sptr0) = 1;
      }

      // process sieve1 (1 mod 4 integers)
      rem(off, LL1, P2);
      if (!IsZero(off))
        rem(off, InvMod(to<ZZ>(4), P2) * to<ZZ>(-off), P2);
      if (off < 0)
        off += P2;

      if (off.SinglePrecision()) {
        conv(off1, off);
        sptr1 = sieve1 + off1;
        if (sptr1 >= sieve1 && sptr1 <= end1)
          (*sptr1) = 1;
      }
    }

    ++pl;
    p = *pl;
    P2 = to<ZZ>(p) * to<ZZ>(p);
  }

  // sieve1 represents squarefree ints that are 1 mod 4
  // sieve0 represetns ints 0 mod 4 with no odd square factors.  Record those
  // for which
  //   x / 4 is 2 or 3 mod 4

  // get first D = 0 mod 4
  PP0 = LL0;
  sptr0 = sieve0;
  while (*sptr0) {
    PP0 += 4;
    ++sptr0;
  }


  // get first D = 1 mod 4
  PP1 = LL1;
  sptr1 = sieve1;
  while (*sptr1) {
    PP1 += 4;
    ++sptr1;
  }


  n = 0;
  while (PP0 <= H || PP1 <= H) {
    if (PP1 < PP0 && PP1 <= H) {
      Dlist[n] = PP1;
      ++n;

      // get next D = 1 mod 4
      do {
        PP1 += 4;
        ++sptr1;
      } while (PP1 <= H && (*sptr1));
    } else {
      off = PP0 >> 2;
      p2l = rem(off, 4);
      if (p2l == 2 || p2l == 3) {
        Dlist[n] = PP0;
        ++n;

      }

      // get next D = 0 mod 4
      do {
        PP0 += 4;
        ++sptr0;
      } while (PP0 <= H && (*sptr0));
    }
  }
}

#endif // AUXILLARYFUNCTIONS_HPP_INCLUDED
