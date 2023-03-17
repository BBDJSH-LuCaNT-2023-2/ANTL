/**
 * @file iq-data.hpp
 * @author Shantha Ramachandran
 *
 *
 */


#ifndef IQ_DATA_H 
#define IQ_DATA_H 

// original global definitions for a local machine, may cause stack overflows

// [MAXP1][MAXP2][MAXP3][MAXP4] = Array of size 40 000
// [MAXPLARGE][MAXE1][MAXE2][MAXE3][MAXE4][MAXE5][MAXE6] = array of size 640 000

#define MAXP 40
#define MAXH 30
#define LASTPRIME 173
#define MAXE1 20
#define MAXE2 10
#define MAXE3 5
#define MAXE4 2
#define MAXE5 2
#define MAXE6 2

#define MAXPHUGE 120
#define MAXPLARGE 80
#define LASTPLARGE 409
#define MAXPSMALL 5

#define MAXP1 10
#define MAXP2 40
#define MAXP3 20
#define MAXP4 5

#define NUMH 10

// global definitions for a local machine
// #define MAXP 5
// #define MAXH 5
// #define LASTPRIME 5
// #define MAXE1 5
// #define MAXE2 5
// #define MAXE3 5
// #define MAXE4 2
// #define MAXE5 2
// #define MAXE6 2
//
// #define MAXPHUGE 5
// // #define MAXPLARGE 80
// #define MAXPLARGE 5
// #define LASTPLARGE 5
// #define MAXPSMALL 5
//
// #define MAXP1 5
// #define MAXP2 5
// #define MAXP3 5
// #define MAXP4 5
//
// #define NUMH 10


#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <boost/multi_array.hpp>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

NTL_CLIENT
using namespace ANTL;

class iq_data {
public:

  ZZ maxD;		// current upper bound on discriminants processed

  long long total0;		// total fields with D = 0 (mod 4)
  long long total1;		// total fields with D = 1 (mod 4)

  long maxp;		// maximum prime required for one class group
  ZZ maxpD;		// its discriminant
  long localmaxp;
  ZZ localmaxpD;

  long maxh0;		// maximum class number D = 0 (mod 4)
  long maxh1;		// maximum class number D = 1 (mod 4)
  ZZ maxh0D;		// the discriminant
  ZZ maxh1D;		// the discriminant


  // *** 1. L(1,X) data ***

  double minL0;		// minimum value of L(1,X) for all D = 0 (mod 4)
  double minL1;		// minimum value of L(1,X) for all D = 1 (mod 8)
  double minL5;		// minimum value of L(1,X) for all D = 5 (mod 8)
  double maxL0;		// maximum value of L(1,X) for all D = 0 (mod 4)
  double maxL1;		// maximum value of L(1,X) for all D = 1 (mod 8)
  double maxL5;		// maximum value of L(1,X) for all D = 5 (mod 8)

  double minLLI0;	// minimum value of LLI for all D = 0 (mod 4)
  double minLLI1;	// minimum value of LLI for all D = 1 (mod 8)
  double minLLI5;	// minimum value of LLI for all D = 5 (mod 8)
  double maxULI0;	// maximum value of ULI for all D = 0 (mod 4)
  double maxULI1;	// maximum value of ULI for all D = 1 (mod 8)
  double maxULI5;	// maximum value of ULI for all D = 5 (mod 8)

  ZZ minLLI0D;
  ZZ minLLI1D;
  ZZ minLLI5D;
  ZZ maxULI0D;
  ZZ maxULI1D;
  ZZ maxULI5D;

  RR sumL0;	        // sum of L(1,X) for all D = 0 (mod 4)
  RR sumL1;	        // sum of L(1,X) for all D = 1 (mod 4)


  // *** 3a. Non-cyclic odd parts of class groups

  long long noncyc0;		// total fields with noncyclic CL and D = 0 (mod 4)
  long long noncyc1;		// total fields with noncyclic CL and D = 1 (mod 4)

  // *** 2. Divisibility by odd primes ***

//   boost::multi_array<long long, 2> pdivs0;
//   boost::multi_array<long long, 2> pdivs1;

  long long pdivs0[MAXP][MAXPSMALL];
  long long pdivs1[MAXP][MAXPSMALL];

  // *** 3b. p-rank probablilities ***

  // prank data
//   boost::multi_array<long long, 7> prank0;
//   boost::multi_array<long long, 7> prank1;
  long long prank0[MAXPLARGE][MAXE1][MAXE2][MAXE3][MAXE4][MAXE5][MAXE6];	// prank data
  long long prank1[MAXPLARGE][MAXE1][MAXE2][MAXE3][MAXE4][MAXE5][MAXE6];	// prank data

  // *** 4. First occurrences of p-sylow groups ***

//   boost::multi_array<ZZ, 7> first_prank0;
//   boost::multi_array<ZZ, 7> first_prank1;
  ZZ first_prank0[MAXPLARGE][MAXE1][MAXE2][MAXE3][MAXE4][MAXE5][MAXE6];	// prank data
  ZZ first_prank1[MAXPLARGE][MAXE1][MAXE2][MAXE3][MAXE4][MAXE5][MAXE6];	// prank data

//   boost::multi_array<long long, 4> two_noncyc0;
//   boost::multi_array<long long, 4> two_noncyc1;

  long long two_noncyc0[MAXP1][MAXP2][MAXP3][MAXP4];	// multiply noncyclic data
  long long two_noncyc1[MAXP1][MAXP2][MAXP3][MAXP4];	// multiply noncyclic data

//   boost::multi_array<ZZ, 4> first_two_noncyc0;
//   boost::multi_array<ZZ, 4> first_two_noncyc1;
  ZZ first_two_noncyc0[MAXP1][MAXP2][MAXP3][MAXP4];	// multiply noncyclic data
  ZZ first_two_noncyc1[MAXP1][MAXP2][MAXP3][MAXP4];	// multiply noncyclic data

  // *** 5. Number of generators ***
  double maxlograt;	// maximum ratio of p/log(Delta)
  double maxloglograt;	// maximum ratio of p/log^2(Delta)
  ZZ maxlogratD;	// the discriminant
  ZZ maxloglogratD;	// the discriminant

  RR sumlograt;	        // sums of ratio of p/log(Delta)
  RR sumloglograt;	// sums of ratio of p/log^2(Delta)

  long long psplit0[MAXPLARGE];      // number of fields needing maxp = p
  long long psplit1[MAXPLARGE];      // number of fields needing maxp = p
  ZZ first_psplit0[MAXPLARGE];  // first field needing maxp = p
  ZZ first_psplit1[MAXPLARGE];  // first field needing maxp = p

  long long kpsplit0[MAXPLARGE];
  long long kpsplit1[MAXPLARGE];
  ZZ first_kpsplit0[MAXPLARGE];
  ZZ first_kpsplit1[MAXPLARGE];


  // *** 6. abc conjecture ***
  long smallh[NUMH];
  ZZ smallhD[NUMH];  


  // *** 7. Odd parts of class numbers *** 
  long maxhodd0;	// maximum odd part of class number D = 0 (mod 4)
  long maxhodd1;	// maximum odd part of class number D = 1 (mod 4)
  ZZ maxhodd0D;	        // the discriminant
  ZZ maxhodd1D;	        // the discriminant

  long ttime;		// total runtime (100th seconds)


  iq_data();
  ~iq_data() {};

  void read_file(char *name, char *pfix);
  void write_file(char *name, char *pfix);

  void write_file_readable(char *name, char *pfix);
//   void write_zz_readable(ZZ &a, FILE *fp);
//   void write_rr_readable(RR &a, FILE *fp);

  void combine(iq_data & newdata);

};


#endif
