/**
 * @file qcl_stats.cpp
 * @author Shantha Ramachandran, Mike Jacobson
 *
 *
 *
 */

#define FILE_INDEX 0
#define PROCESS 1
#define NO_PROCESS 2
#define FILE_DONE 3
#define SLAVE_DONE 4

// #define LIST_SIZE_QUADRATIC 33554432
#define LIST_SIZE_QUADRATIC 10000

#include "iq-data.hpp"
#include "iq-stats.hpp"
#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <dirent.h>
#include <fstream>
#include <list>
#include <mpi.h>
#include <unistd.h>

#include "AuxillaryFunctions.hpp"
#include "timer.hpp"

NTL_CLIENT
using namespace ANTL;

int process_file(int index, char *fname, char *iprefix, char *oprefix);
void output_final(iq_data data, long machines);
int p_rank(long CL[], int rank, long p);

double nu(long p);
double nu(long p, long alpha);

void create_output_files();
int output_merged(int numfiles, iq_data data);
int output_interval(iq_data data, ZZ thresh);
int output_graph_interval(iq_data data, ZZ gthresh);
void close_output_files();

// **********************************************************
// *** Main program to compute statistics on class groups ***
// **********************************************************
int main(int argc, char **argv) {

  std::cout << "HARVST: HELOOOOOOOOOO" << std::endl;
  bool DBG_HARVST = true;
  if(DBG_HARVST) {std::cout << "HARVST: BEGIN" << std::endl;};
  // parameters
  int myrank;

  // initialize mpi
  if(DBG_HARVST) {std::cout << "HARVST: INIT" << std::endl;};
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // **********************
  // *** master process ***
  // **********************
  if (myrank == 0) {
    std::cout << "MASTER: BEGIN" << std::endl;

    if(DBG_HARVST) {std::cout << "MASTER: INITIALIZE VARIABLES" << std::endl;};
    // master parameters
    int done, machines, i, j, k, index;
    long L, H, bound, current, total;
    ZZ IS, total_time;
    char hname[50], fname[50], back[50];
    char iprefix[50], oprefix[50];
    ifstream infile;
    ofstream outfile;
    MPI_Status status;
    iq_data combined_data;
    iq_data new_data;
    ZZ threshold, gthreshold, inc, ginc;
    ZZ compL, compH;
    timer t1;

    // parameters for interval computation
    ZZ Dlist[LIST_SIZE_QUADRATIC];
//     quadratic_order<long> QOl;
//     quadratic_order<long long> QO;
    vec_ZZ Cl;
//     long Dl, oldDl, rank, n;
//     long long D, oldD;
    long n;
    ZZ Dl, oldDl, rank;
    ZZ D, oldD;
    char zipper[50], syscall[50];

    // get number of available machines
    MPI_Comm_size(MPI_COMM_WORLD, &machines);
    machines--;

    // check arguments
    if (argc != 6) {
      cerr << "Usage: iq_stats L H IS iprefix oprefix" << endl;
      cerr << "L - lower index for interval" << endl;
      cerr << "H - upper index for interval" << endl;
      cerr << "IS - interval size" << endl;
      cerr << "iprefix - prefix for nodes' local data files" << endl;
      cerr << "oprefix - prefix for output data files" << endl;
      exit(1);
    }

    if(DBG_HARVST) {std::cout << "MASTER: GET ARGUEMENTS" << std::endl;};
    // get arguments
    L = atol(argv[1]);
    H = atol(argv[2]);
    conv(IS, argv[3]);
    sprintf(iprefix, "%s", argv[4]);
    sprintf(oprefix, "%s", argv[5]);

    if(DBG_HARVST) {std::cout << "MASTER: PRINT OUT MACHINE INFO" << std::endl;};
    // print out machine information
    gethostname(hname, 50);
    cout << "Master: " << hname << endl;
    cout << "Computing statistics for all class numbers for fundamental D "
         << "in interval " << to<ZZ>(L) * IS << " to " << to<ZZ>(H) * IS
         << endl;
    cout << "Using " << machines << " processes" << endl;

    // create interval thresholds
    inc = 10000000;
    ginc = 10000000;
    threshold = inc;
    gthreshold = ginc;

    // create checkpoint table
    int finished[H - L];
    for (i = 0; i < H - L; i++) {
      finished[i] = 0;
    }

    if(DBG_HARVST) {std::cout << "MASTER: LOAD CHECKPOINT TABLE (IF IT EXISTS)" << std::endl;};
    // load checkpoint table (if exists)
    sprintf(fname, "stats-check-%ld-%s", L, argv[3]);
    cout << "Checkpoint file = " << fname << endl;
    infile.open(fname);

    if (infile.fail()) {
      // clear values
      clear(total_time);
      total = 0;
      current = 0;
      infile.clear();
    } else {
      // read in values
      infile >> total_time;
      infile >> total;
      infile >> current;
      infile >> threshold;
      infile >> inc;
      infile >> gthreshold;
      for (i = 0; i < H - L; i++) {
        infile >> finished[i];
        if (finished[i] == 1)
          finished[i] = 0;
      }

      infile.close();
    }

    if (current > 0) {
      // create combined output file
      combined_data.read_file("/current.dat", oprefix);
    } else {
      // create latex output files
      create_output_files();
    }

    if(DBG_HARVST) {std::cout << "MASTER: INITIALIZE PRIME LIST" << std::endl;};
    // initialize prime list
//     get_Dlist_imag(to<ZZ>(L), to<ZZ>(H * IS), Dlist, n, 1);
    get_Dlist_real(to<ZZ>(L), to<ZZ>(H * IS), Dlist, n, 1);

    if(DBG_HARVST) {std::cout << "MASTER: GATHER STATS" << std::endl;};
    // gather stats until all processes report finished
    done = 0;
    while (done != machines) {

      // receive message from slave
      MPI_Recv(&index, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
               &status);

      if (status.MPI_TAG == FILE_INDEX) {

        // check if file has been assigned
        if (finished[index] == 0) {
          // if file has not been assigned, send a message to process it
          MPI_Send(&myrank, 1, MPI_INT, status.MPI_SOURCE, PROCESS,
                   MPI_COMM_WORLD);

          // assign file in table
          finished[index] = 1;

        } else {
          // if file has been assigned, send a message to leave it
          MPI_Send(&myrank, 1, MPI_INT, status.MPI_SOURCE, NO_PROCESS,
                   MPI_COMM_WORLD);

        } // end if else

      } else if (status.MPI_TAG == FILE_DONE) {

        // mark file as being done
        finished[index] = 2;
        total++;

        // update checkpoint file
        outfile.clear();
        sprintf(fname, "stats-check-%ld-%s", L, argv[3]);
        outfile.open(fname);
        if (outfile.fail()) {
          cerr << " Could not open check file " << fname << endl;
        }

        outfile << total_time << " " << total << " " << current << " ";
        outfile << threshold << " " << inc << " " << gthreshold << " ";
        for (i = 0; i < H - L; i++) {
          outfile << finished[i] << " ";
        }
        outfile << endl;
        outfile.close();

      } else if (status.MPI_TAG == SLAVE_DONE) {

        // one more slave has finished
        done++;

      } // end if else

    } // end while

    cout << "Combining files ..." << endl;

    // check to make sure all intervals have been processed
    long oldcurrent = current;
    for (i = current; i < H; i++) {
      if (finished[i] != 2) {
        cout << "Missing interval " << i << endl;

        // compute and process this interval
        compL = to<ZZ>(IS * i + 1);
        compH = to<ZZ>(IS * (i + 1));
        cout << "Computing interval from " << compL << " to " << compH << endl;

        // create file
        sprintf(fname, "%s/iq-data.%d", iprefix, i);
        outfile.open(fname);

        t1.start_timer();

        // get Dlist
//         get_Dlist_imag(compL, compH, Dlist, n);
//         get_Dlist_real(compL, compH, Dlist, n, 0);

        if (compL == 1) {
          get_Dlist_real(to_ZZ(2), compH, Dlist, n, 0);
        }
        else {
          get_Dlist_real(compL, compH, Dlist, n, 0);
        }

        outfile << compL << " " << compH << " " << n << endl;

        // set global # of terms and prec for L-function approximations for
        // this interval
//         QO.set_Lfunction_global(to<long long>(compH));

        oldD = to<ZZ>(compL);
        for (j = 0; j < n; ++j) {
          D = to<ZZ>(Dlist[j]);
          std::cout << "Current discriminant is " << D << std::endl;

//           QO.assign(D);
          QuadraticOrder<ZZ> QO{ZZ(D)};
//           Cl = QO.class_group(CLASS_GROUP_BSGS);
          pair<double, vector<long>> regulator_and_class_group = get_regulator_and_class_group(QO);
          double regulator = regulator_and_class_group.first;
          vector<long> class_group = regulator_and_class_group.second;

//           rank = QO.get_rank();
          rank = class_group.size();

          outfile << (D - oldD) << " " << flush;
//           outfile << QO.get_nump() << " " << flush;
//           outfile << QO.get_pmax() << " " << flush;
          outfile << 1 << " " << 1 << " " << flush;

          outfile << regulator << " " << flush;
          outfile << rank << flush;
          for (k = 0; k < rank; ++k)
            outfile << " " << class_group[k] << flush;
          outfile << endl;
//           oldD = ANTL::abs(D);
          oldD = D;
        }

        // get time for interval
        t1.stop_timer();
        outfile << t1.user_time() << endl;
        outfile.close();
        std::cout << "TEST 1" << std::endl;
        // zip file
        sprintf(zipper, "gzip %s/iq-data.%d", iprefix, i);
        system(zipper);
        std::cout << "TEST 2" << std::endl;
        // process interval
        sprintf(fname, "%s/iq-data.%d.gz", iprefix, i);
        process_file(i, fname, iprefix, oprefix);
        std::cout << "TEST 3" << std::endl;

        // mark file as being done
        finished[i] = 2;
        total++;

        std::cout << "TEST 4" << std::endl;
        // update checkpoint file
        outfile.clear();
        sprintf(fname, "stats-check-%ld-%s", L, argv[3]);
        outfile.open(fname);
        if (outfile.fail()) {
          cerr << " Could not open check file " << fname << endl;
        }

        std::cout << "TEST 5" << std::endl;
        outfile << total_time << " " << total << " " << oldcurrent << " ";
        outfile << threshold << " " << inc << " " << gthreshold << " ";
        for (j = 0; j < H - L; j++) {
          outfile << finished[j] << " ";
        }
        outfile << endl;
        outfile.close();
      }

      // unzip the interval data
      sprintf(syscall, "gunzip %s/data%d.dat.gz", oprefix, i);
      system(syscall);

      // read the interval data
      sprintf(fname, "/data%d.dat", i);
      new_data.read_file(fname, oprefix);

      // zip the interval data back up
      sprintf(syscall, "gzip %s/data%d.dat", oprefix, i);
      system(syscall);

      // combine the new data
      combined_data.combine(new_data);

      // output the threshold data
      if (combined_data.maxD > threshold - 100) {

        if (output_interval(combined_data, threshold) == 1) {
          cerr << "Could not output interval " << threshold << "-"
               << threshold + inc << endl;
          exit(1);
        }
        threshold += inc;
        if (threshold == 10 * inc) {
          inc *= 10;
        }

        for (j = 0; j < NUMH; j++) {
          combined_data.smallh[j] = 0;
          combined_data.smallhD[j] = 0;
        }
      }

      // output the graph threshold data
      if (combined_data.maxD > gthreshold - 100) {
        if (output_graph_interval(combined_data, gthreshold) == 1) {
          cerr << "Could not output interval " << gthreshold << "-"
               << gthreshold + ginc << endl;
          exit(1);
        }
        gthreshold += ginc;
      }

      current++;
      std::cout << "GOT TO HERE" << std::endl;
    }

    // write the total combined file
    combined_data.write_file("/current.dat", oprefix);

    // back up the combined file
    sprintf(back, "/back%ld.dat", current);
    combined_data.write_file(back, oprefix);

    // update checkpoint file
    outfile.clear();
    sprintf(fname, "stats-check-%ld-%s", L, argv[3]);
    outfile.open(fname);
    outfile << total_time << " " << total << " " << current << " ";
    outfile << threshold << " " << inc << " " << gthreshold << " ";
    for (i = 0; i < H - L; i++) {
      outfile << finished[i] << " ";
    }
    outfile << endl;
    outfile.close();

    // merge some results
    output_merged(H - L, combined_data);

    // output some final results
    output_final(combined_data, machines);

    // close all the files
    close_output_files();

  } // end master process
  // *********************
  // *** slave process ***
  // *********************
  else {

    if(DBG_HARVST) {std::cout << "SLAVE: BEGIN" << std::endl;};
    if(DBG_HARVST) {std::cout << "SLAVE: INITIALIZE VARIABLES" << std::endl;};
    // slave parameters
    char hname[50], fname[50], dname[50];
    char iprefix[50], oprefix[50];
    char fsub[50], buffer[50];
    char *place;
    DIR *dp;
    struct dirent *dirp;
    int index;
    MPI_Status status;
    long long L, H;
    long bound;

    L = atol(argv[1]);
    H = atol(argv[2]);
    sprintf(iprefix, "%s", argv[4]);
    sprintf(oprefix, "%s", argv[5]);

    // print out process information
    gethostname(hname, 50);
    cout << "Process " << myrank << " (" << hname << "): Initializing" << endl;

    if(DBG_HARVST) {std::cout << "SLAVE: FIND OUT HOST FILES" << std::endl;};
    // find out what files exist on the host
    sprintf(fsub, "iq-data.");
    sprintf(dname, iprefix);
    if (!(dp = opendir(dname))) {
      cout << myrank << " Can't open directory: " << dname << endl;

    } else {

    if(DBG_HARVST) {std::cout << "SLAVE: FOR EACH HOST FILE" << std::endl;};
      // for each file on host
      while ((dirp = readdir(dp)) != 0) {

        // check if file is correct format
        sprintf(fname, dirp->d_name);
        place = strstr(fname, fsub);
        if (place) {

          // get index of file
          place = strstr(fname, ".") + 1;
          index = atoi(place);

          if (index < H && index >= L) {

            // send message to master
            MPI_Send(&index, 1, MPI_INT, 0, FILE_INDEX, MPI_COMM_WORLD);

            // receive ack from master
            MPI_Recv(&buffer, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD,
                     &status);

            // only process the file if the master tells us to
            if (status.MPI_TAG == PROCESS) {

              // process file
              strcpy(buffer, dname);
              strcat(buffer, "/");
              strcat(buffer, fname);

              cout << "Process " << myrank
                   << " computing statistics on interval " << index << endl;
              if (process_file(index, buffer, dname, oprefix)) {
                cerr << "Error processing file" << endl;
                // exit(1);
              } else {
                // send done index to master
                MPI_Send(&index, 1, MPI_INT, 0, FILE_DONE, MPI_COMM_WORLD);
              }

            } // end if

          } // end if

        } // end if

      } // end while

      // close directory
      closedir(dp);
    }

    // send finished message to master
    MPI_Send(&myrank, 1, MPI_INT, 0, SLAVE_DONE, MPI_COMM_WORLD);

    cout << "Process " << myrank << ": Finished" << endl;
  }/*

  } // end slave
  */

  // finalize mpi
  MPI_Finalize();
  if (myrank == 0) {
    cout << "Finished statistics gathering" << endl;
  }

  return 0;
} // end main

// ****************************************************
// *** function to process a class number data file ***
// ****************************************************
int process_file(int index, char *fname, char *iprefix, char *oprefix) {
  std::cout << "process_file: begin!" << std::endl;
  // parameters
  ifstream infile;
  ofstream tempfile, errorfile;
  ZZ L, H, offset;
  int prank;
  long n, rank, num_p, max_p, i, j, k, l;
  char zipper[50], unzipped[50], datatemp[50], syscall[50];
  ZZ D, r;
  long h, hodd, p;
  RR L1X, LD, temp, maxrat;
  PrimeSeq sp;
  iq_data data;
  timer time;
  int e1, e2, e3, e4, e5, e6;
  int ex[6];
  int D4;
  long temp2;
  RR c, LLI, ULI;
  long multp[MAXPLARGE];
  long numnoncyc;
  int numsmallh;
  int pidx[LASTPLARGE];

  list<pair<long, ZZ>> smallh;
  list<pair<long, ZZ>>::iterator iter;

  // start timer
  time.start_timer();

  // unzip file
  sprintf(zipper, "gunzip %s", fname);
  if (system(zipper) != 0)
    perror("Could not unzip iq file");

  // find unzipped file
  sprintf(unzipped, "%s/iq-data.%d", iprefix, index);

  // open file
  infile.open(unzipped);
  if (infile.fail()) {
    cerr << "ERROR: File " << unzipped << " missing" << endl << flush;
    return 1;
  }

  // open data file for L(1,X)
  sprintf(datatemp, "%s%d.tex", datafile, index);
  tempfile.open(datatemp, ofstream::out);
  if (tempfile.fail()) {
    cerr << "Error opening " << datatemp << endl;
    return 1;
  }

  // initializations
  std::cout << "process_file: initializations" << std::endl;
  numsmallh = 0;
  sp.reset(3);
  for (i = 0; i < MAXPLARGE; i++) {
    p = sp.next();
    pidx[p] = i;
  }

  // read header values from file
  std::cout << "process_file: read header values" << std::endl;
  infile >> L;
  infile >> H;
  infile >> n;
  std::cout << "are header values" << L << H << n << std::endl;

  // read data from file
  std::cout << "process_file: read data from file" << std::endl;
  D = L;
  for (i = 0; i < n; i++) {
    std::cout << "iterating!" << std::endl;
    // calculate new discriminant D
    std::cout << "calculate new discriminant D" << std::endl;
    infile >> offset;
    D += offset;
    data.maxD = D;
    std::cout << "new discriminant D is " << D << std::endl;

    // get number of prime generators
    std::cout << "get number of prime generators" << std::endl;
    infile >> num_p;

    // get max prime
    std::cout << "get max prime" << std::endl;
    infile >> max_p;
    if (max_p > data.maxp) {
      data.maxp = max_p;
      data.maxpD = D;
      data.localmaxp = max_p;
      data.localmaxpD = D;
    }

    // get regulator
    double regulator;
    infile >> regulator;
    std::cout << "regulator is " << regulator << std::endl;

    // get rank of field and CL
    std::cout << "get rank of field and CL" << std::endl;
    infile >> rank;
    std::cout << "rank is " << rank << std::endl;
    long Cl[rank];
    h = 1;
    for (j = 0; j < rank; j++) {
      infile >> Cl[j];
      h *= Cl[j];
    }

    // find small h
    std::cout << "find small h" << std::endl;
    pair<long, ZZ> maxsmall;
    if (smallh.size() < 10) {
      pair<long, ZZ> newpair(h, D);
      smallh.push_back(newpair);
      smallh.sort();
    } else {
      maxsmall = smallh.back();
      if (maxsmall.first > h) {
        pair<long, ZZ> newpair(h, D);
        smallh.insert(smallh.end()--, newpair);
        smallh.pop_back();
        smallh.sort();
      }
    }

    // calculate hodd
    std::cout << "calculate hodd" << std::endl;
    hodd = h;
    while (!hodd & 1)
      hodd <<= 1;

    // calculate odd part of CL
    std::cout << "calculate odd part of CL" << std::endl;
    long oddCl[rank];
    for (j = 0; j < rank; j++) {
      oddCl[j] = Cl[j];
      if ((oddCl[j] & 1) == 0)
        oddCl[j] >>= 1;
    }

    // calculate L(1,X)
    std::cout << "calculate L(1,X)" << std::endl;
    L1X = to<RR>(h) * ComputePi_RR() / sqrt(to<RR>(D));

    // compute max ideal stats
    std::cout << "compute max ideal stats" << std::endl;
    if (D > 4) {
      temp = log(D);
      maxrat = to<RR>(max_p) / temp;
      data.sumlograt += maxrat;
      if (maxrat > data.maxlograt) {
        data.maxlograt = to<double>(maxrat);
        data.maxlogratD = D;
      }
      temp *= temp;
      maxrat = to<RR>(max_p) / temp;
      data.sumloglograt += maxrat;
      if (maxrat > data.maxloglograt) {
        data.maxloglograt = to<double>(maxrat);
        data.maxloglogratD = D;
      }
    }

    // find value (mod 4)
    rem(r, D, to<ZZ>(8));
    if (r == 0 || r == 4)
      D4 = 1;
    else
      D4 = 0;

    // look for max p splits
    std::cout << "look for max p splits" << std::endl;
    if (max_p > LASTPLARGE || num_p > MAXPLARGE) {
      cout << "Prime Ideal: D = " << D << " max p = " << max_p
           << " num p = " << num_p << endl;
    } else {

      if (D4) {
        if (data.psplit0[pidx[max_p]] == 0)
          data.first_psplit0[pidx[max_p]] = D;
        data.psplit0[pidx[max_p]]++;

        if (data.kpsplit0[num_p] == 0)
          data.first_kpsplit0[num_p] = D;
        data.kpsplit0[num_p]++;

      } else {
        if (data.psplit1[pidx[max_p]] == 0)
          data.first_psplit1[pidx[max_p]] = D;
        data.psplit1[pidx[max_p]]++;

        if (data.kpsplit1[num_p] == 0)
          data.first_kpsplit1[num_p] = D;
        data.kpsplit1[num_p]++;
      }
    }

    // compute prank data
    std::cout << "compute prank data" << std::endl;
    numnoncyc = 0;
    for (j = 0; j < MAXPLARGE; j++)
      multp[j] = 0;

    sp.reset(2);
    long second = oddCl[1];
    int donep = 0;
    if (rank == 1 || oddCl[1] == 1)
      donep = 1;

    j = 0;
    while (!donep || j < MAXP) {
      std::cout << "im trapped here" << std::endl;
      p = sp.next();
      prank = p_rank(oddCl, rank, p);

      if (prank > 6) {
        cout << "P-rank: D = " << D << " " << p << "-rank = " << prank << endl;
      }

      if (prank > 1) {

        if (j >= MAXPLARGE) {
          cout << "P-rank: D = " << D << " " << p << "-rank = " << prank
               << endl;
        } else {

          // save prank
          multp[numnoncyc] = j;
          numnoncyc++;

          // compute exponents for p-rank
          ex[0] = ex[1] = ex[2] = ex[3] = ex[4] = ex[5] = 0;
          l = 0;

          k = 0;
          while ((k < rank) && (l < 6)) {
            ex[l] = 0;
            temp2 = oddCl[k];
            while (!(temp2 % p)) {
              ex[l]++;
              temp2 /= p;
            }
            l++;
            k++;
          }
          e1 = ex[0];
          e2 = ex[1];
          e3 = ex[2];
          e4 = ex[3];
          e5 = ex[4];
          e6 = ex[5];

          if (e1 >= MAXE1 || e2 >= MAXE2 || e3 >= MAXE3 || e4 >= MAXE4 ||
              e5 >= MAXE5 || e6 >= MAXE6) {
            cout << "P-rank: D = " << D << ": " << p << "-rank = " << prank
                 << endl;
            cout << "Exponents = " << e1 << ' ' << e2 << ' ' << e3 << ' ' << e4
                 << ' ' << e5 << ' ' << e6 << endl;
          } else {

            // save prank data
            if (D4) {
              if (data.prank0[j][e1][e2][e3][e4][e5][e6] == 0)
                data.first_prank0[j][e1][e2][e3][e4][e5][e6] = D;
              data.prank0[j][e1][e2][e3][e4][e5][e6]++;
            } else {
              if (data.prank1[j][e1][e2][e3][e4][e5][e6] == 0)
                data.first_prank1[j][e1][e2][e3][e4][e5][e6] = D;
              data.prank1[j][e1][e2][e3][e4][e5][e6]++;
            }
          }
        }
        // second /= p;
        while (second % p == 0)
          second /= p;
        if (second == 1) {
          donep = 1;
        }
      }

      if (j < MAXP) {
        // save divisor of h data
        if (h % p == 0) {
          if (D4)
            data.pdivs0[j][0]++;
          else {
            data.pdivs1[j][0]++;
          }
        }

        // check if p^i | h
        long pp;
        pp = p;
        for (k = 1; k < MAXPSMALL; k++) {
          pp *= p;
          if (h % pp == 0) {
            if (D4)
              data.pdivs0[j][k]++;
            else
              data.pdivs1[j][k]++;
          }
        }
      }
      j++;

      if(j > 1000) break;
    }

    // doubly and trebly noncyclic
    std::cout << "doubly and trebly noncyclic" << std::endl;
    if (numnoncyc > 1) {

      if (numnoncyc > 6) {
        cout << "Multiply Non-cyclic: D = " << D << " noncyclic > 6" << endl;
      }

      ex[0] = ex[1] = ex[2] = ex[3] = ex[4] = ex[5] = 0;
      for (j = 0; j < numnoncyc; j++) {
        ex[j] = multp[j];
      }
      e1 = ex[0];
      e2 = ex[1];
      e3 = ex[2];
      e4 = ex[3];
      if (ex[4] > 0) {
        cout << "Multiply Non-cyclic: D = " << D << "more noncyclic " << endl;
        cout << "Exponents: " << e1 << ' ' << e2 << ' ' << e3 << ' ' << e4
             << ' ' << ex[4] << ' ' << ex[5] << endl;
      }
      if (e1 >= MAXP1 || e2 >= MAXP2 || e3 >= MAXP3 || e4 >= MAXP4) {
        cout << "Multiply Non-cyclic: D = " << D << endl;
        cout << "Exponents: " << e1 << ' ' << e2 << ' ' << e3 << ' ' << e4
             << ' ' << ex[4] << ' ' << ex[5] << endl;
      } else {

        if (D4) {
          if (data.two_noncyc0[e1][e2][e3][e4] == 0)
            data.first_two_noncyc0[e1][e2][e3][e4] = D;
          data.two_noncyc0[e1][e2][e3][e4]++;
        } else {
          if (data.two_noncyc1[e1][e2][e3][e4] == 0)
            data.first_two_noncyc1[e1][e2][e3][e4] = D;
          data.two_noncyc1[e1][e2][e3][e4]++;
        }
      }
    }

    // check if oddCL is noncyclic
    std::cout << "check if oddCL is noncyclic" << std::endl;
    if (numnoncyc > 0 && multp[0] != 0) {
      if (D4)
        data.noncyc0++;
      else
        data.noncyc1++;
    }

    if (D4) {
      // fields with D = 0 (mod 4)

      // update total fields with D = 0 (mod 4)
      data.total0++;

      // update min and max values of L(1,X)
      if (L1X < data.minL0) {
        data.minL0 = to<double>(L1X);
        tempfile << minL0_c << ' ' << D << ' ' << L1X << endl;
      }
      if (L1X > data.maxL0) {
        data.maxL0 = to<double>(L1X);
        tempfile << maxL0_c << ' ' << D << ' ' << L1X << endl;
      }
      data.sumL0 += L1X;

      // calculate LLI and ULI
      c = 8 * ANTL::exp(gamma) / (ComputePi_RR() * ComputePi_RR());
      LLI = L1X * c * log(log(to_RR(D)));
      c = ANTL::exp(gamma);
      ULI = L1X / (c * log(log(to_RR(D))));

      // update LLI, ULI
      if (LLI < data.minLLI0) {
        if (D != 4 && D != 232) {
          data.minLLI0 = to<double>(LLI);
          data.minLLI0D = D;
          tempfile << minLLI0_c << ' ' << D << ' ' << L1X << endl;
        }
      }
      if (ULI > data.maxULI0) {
        if (D != 4 && D != 8 && D != 20) {
          data.maxULI0 = to<double>(ULI);
          data.maxULI0D = D;
          tempfile << maxULI0_c << ' ' << D << ' ' << L1X << endl;
        }
      }

      // update max class number
      if (h > data.maxh0) {
        data.maxh0 = h;
        data.maxh0D = D;
      }

      // odd h data
      if (hodd > data.maxhodd0) {
        data.maxhodd0 = hodd;
        data.maxhodd0D = D;
      }

    } else {

      // update total fields with D = 1 (mod 4)
      data.total1++;

      // update L(1,X)
      data.sumL1 += L1X;

      // calculate LLI and ULI
      c = 12 * ANTL::exp(gamma) / (ComputePi_RR() * ComputePi_RR());
      LLI = L1X * c * log(log(to_RR(D)));
      c = 2 * ANTL::exp(gamma);
      ULI = L1X / (c * log(log(to_RR(D))));

      // update max class number
      if (h > data.maxh1) {
        data.maxh1 = h;
        data.maxh1D = D;
      }

      // odd h data
      if (hodd > data.maxhodd1) {
        data.maxhodd1 = hodd;
        data.maxhodd1D = D;
      }

      if (r == 7) {
        // fields with D = 1 (mod 8)

        // update min and max values of L(1,X)
        if (L1X < data.minL1) {
          data.minL1 = to<double>(L1X);
          tempfile << minL1_c << ' ' << D << ' ' << L1X << endl;
        }
        if (L1X > data.maxL1) {
          data.maxL1 = to<double>(L1X);
          tempfile << maxL1_c << ' ' << D << ' ' << L1X << endl;
        }

        // update ULI, LLI
        if (LLI < data.minLLI1) {
          if (D != 7) {
            data.minLLI1 = to<double>(LLI);
            data.minLLI1D = D;
            tempfile << minLLI1_c << ' ' << D << ' ' << L1X << endl;
          }
        }
        if (ULI > data.maxULI1) {
          data.maxULI1 = to<double>(ULI);
          data.maxULI1D = D;
          tempfile << maxULI1_c << ' ' << D << ' ' << L1X << endl;
        }

      } else {
        // fields with D = 5 (mod 8)

        // update min and max values of L(1,X)
        if (L1X < data.minL5) {
          data.minL5 = to<double>(L1X);
          tempfile << minL5_c << ' ' << D << ' ' << L1X << endl;
        }
        if (L1X > data.maxL5) {
          data.maxL5 = to<double>(L1X);
          tempfile << maxL5_c << ' ' << D << ' ' << L1X << endl;
        }

        // update ULI, LLI
        if (LLI < data.minLLI5) {
          if (D != 3 && D != 163) {
            data.minLLI5 = to<double>(LLI);
            data.minLLI5D = D;
            tempfile << minLLI5_c << ' ' << D << ' ' << L1X << endl;
          }
        }
        if (ULI > data.maxULI5) {
          if (D != 3 && D != 11) {
            data.maxULI5 = to<double>(ULI);
            data.maxULI5D = D;
            tempfile << maxULI5_c << ' ' << D << ' ' << L1X << endl;
          }
        }
      }
    }
  }

  std::cout << "process_file: closing temp file!" << std::endl;
  tempfile.close();

  // update smallh
  std::cout << "process_file: update smallh" << std::endl;
  pair<long, ZZ> newpair;
  long element, tempelement;
  ZZ elementD, tempD;
  for (iter = smallh.begin(); iter != smallh.end(); iter++) {
    newpair = *iter;
    if (newpair.second != 0) {
      element = newpair.first;
      elementD = newpair.second;
      j = 0;
      while (data.smallh[j] < element && data.smallhD[j] != 0)
        j++;
      while (j < 10) {
        tempelement = data.smallh[j];
        tempD = data.smallhD[j];
        data.smallh[j] = element;
        data.smallhD[j] = elementD;
        element = tempelement;
        elementD = tempD;
        j++;
      }
    }
  }

  std::cout << "process_file: stop timer" << std::endl;
  time.stop_timer();
  data.ttime = time.user_time();
  cout << "Time for interval " << index << ": ";
  MyTime(time.user_time());
  cout << endl;

  // close file
  std::cout << "process_file: close file" << std::endl;
  infile.close();

  // write the computed interval data to an output file
  std::cout << "process_file: write the computed interval data to an output file" << std::endl;
  char name[50];
  sprintf(name, "/data%d.dat", index);
  data.write_file(name, oprefix);

  sprintf(syscall, "gzip %s/data%d.dat", oprefix, index);
  if (system(syscall) != 0)
    perror("Could not zip data file");

  sprintf(zipper, "gzip %s", unzipped);
  if (system(zipper) != 0)
    perror("Could not zip iq file");

  std::cout << "process_file: finish!" << std::endl;
  return 0;
}

void output_final(iq_data data, long machines) {

  cout << endl << endl << "*** Final results ***" << endl;

  cout << "Maximum prime required: " << data.maxp << endl;
  cout << "Discriminant: " << data.maxpD << endl << endl;

  cout << "Maximum p/log(Delta): " << data.maxlograt << endl;
  cout << "Discriminant: " << data.maxlogratD << endl << endl;

  cout << "Maximum p/log^2(Delta): " << data.maxloglograt << endl;
  cout << "Discriminant: " << data.maxloglogratD << endl << endl;

  cout << "Maximum h (D = 0 (mod 4)): " << data.maxh0 << endl;
  cout << "Disctiminant: " << data.maxh0D << endl << endl;

  cout << "Maximum h (D = 1 (mod 4)): " << data.maxh1 << endl;
  cout << "Discriminant: " << data.maxh1D << endl << endl;

  cout << "Maximum odd h (D = 0 (mod 4)): " << data.maxhodd0 << endl;
  cout << "Discriminant: " << data.maxhodd0D << endl << endl;

  cout << "Maximum odd h (D = 1 (mod 4)): " << data.maxhodd1 << endl;
  cout << "Discriminant: " << data.maxhodd1D << endl << endl;

  cout << "Total time: ";
  MyTime(data.ttime);
  cout << endl;
  cout << "Real time: ";
  MyTime(data.ttime / machines);
  cout << endl << "*********************" << endl << endl;
}

int p_rank(long CL[], int rank, long p) {

  int prank = 0;
  int i;

  for (i = 0; i < rank; i++) {
    if (CL[i] % p == 0)
      prank++;
  }

  return prank;
}

// function to create latex output files
void create_output_files() {

  char fname[100];
  FILE *fp;

  cout << "Creating output files ..." << endl;

  // **********************
  // *** 1. L(1,x) data ***
  // **********************

  // successive L(1) minima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, minL0file);
  std::cout << "create_o_f: 1 fname is: " << fname << std::endl;
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  std::cout << "create_o_f: 2" << std::endl;

  fprintf(fp, "Successive L(1) minima:\n");
  fprintf(fp, "D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive L(1) minima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minL1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive L(1) minima:\n");
  fprintf(fp, "D = 1 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive L(1) minima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minL5file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive L(1) minima:\n");
  fprintf(fp, "D = 5 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive L(1) maxima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, maxL0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive L(1) maxima:\n");
  fprintf(fp, "D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive L(1) maxima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxL1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive L(1) maxima:\n");
  fprintf(fp, "D = 1 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive L(1) maxima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxL5file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive L(1) maxima:\n");
  fprintf(fp, "D = 5 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive LLI minima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, minLLI0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive LLI minima:\n");
  fprintf(fp, "D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive LLI minima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minLLI1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive LLI minima:\n");
  fprintf(fp, "D = 1 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive LLI minima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minLLI5file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive LLI minima:\n");
  fprintf(fp, "D = 5 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$LLI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive ULI maxima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, maxULI0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive ULI maxima:\n");
  fprintf(fp, "D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive ULI maxima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxULI1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive ULI maxima:\n");
  fprintf(fp, "D = 1 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // successive ULI maxima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxULI5file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Successive ULI maxima:\n");
  fprintf(fp, "D = 5 (mod 8)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L(1,\\chi)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$L_D(1)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ULI_D$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // local LLI minima
  strcpy(fname, fprefix);
  strcat(fname, localLLIgraph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Local LLI minima:\n");
  fclose(fp);

  // local ULI maxima
  strcpy(fname, fprefix);
  strcat(fname, localULIgraph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Local ULI maxima:\n");
  fclose(fp);

  // file of average L function
  strcpy(fname, fprefix);
  strcat(fname, aveLfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Average L\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{overall} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta \\equiv \\Mod{0} {4}$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta \\equiv \\Mod{1} {4}$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // ***********************************
  // *** 2. divisors of h statistics ***
  // ***********************************

  // file of divisors of h statistics
  strcpy(fname, fprefix);
  strcat(fname, pdivsfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{2-8}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$19$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h statistics (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivs0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{2-8}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$19$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h statistics (D = 1 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivs1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{2-8}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$19$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h percentages
  strcpy(fname, fprefix);
  strcat(fname, pdivsperfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (percentages)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{2-8}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$19$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h percentages (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivsper0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp,
          "Class number divisor statistics (percentages) D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{2-8}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$19$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h percentages (D = 1 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivsper1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp,
          "Class number divisor statistics  (percentages) D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{2-8}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$19$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivsratfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics  (ratios)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_3(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_5(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_7(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h ratios, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivsrat0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics  (ratios), D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_3(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_5(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_7(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisors of h ratios, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivsrat1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics  (ratios), D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_3(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_5(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_7(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of square divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivssqrfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (squares)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{3^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{5^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{7^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19^2}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of square divisors of h ratios, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqr0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (squares), D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{3^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{5^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{7^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19^2}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of square divisors of h ratios, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqr1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (squares), D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{3^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{5^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{7^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17^2}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19^2}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of cube divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivscubfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (cubes)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{3^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{5^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{7^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19^3}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of cube divisors of h ratios, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivscub0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (cubes), D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{3^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{5^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{7^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19^3}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of cube divisors of h ratios, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivscub1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Class number divisor statistics (cubes), D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{3^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{5^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{7^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{11^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{13^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{17^3}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_{19^3}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of divisor graph data for p = 3
  strcpy(fname, fprefix);
  strcat(fname, pdivs3graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  fprintf(fp, "# Class number divisor statistics for p = 3\n\n");
  fprintf(fp, "#x\t#p_p(x)\n");
  fclose(fp);

  // file of divisor graph data for p = 5
  strcpy(fname, fprefix);
  strcat(fname, pdivs5graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  fprintf(fp, "# Class number divisor statistics for p = 5\n\n");
  fprintf(fp, "#x\t#p_p(x)\n");
  fclose(fp);

  // file of divisor graph data for p = 5
  strcpy(fname, fprefix);
  strcat(fname, pdivs7graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  fprintf(fp, "# Class number divisor statistics for p = 7\n\n");
  fprintf(fp, "#x\t#p_p(x)\n");
  fclose(fp);

  // *********************************************
  // *** 3a. Noncyclic odd part of class group ***
  // *********************************************

  // file of maximum class numbers
  strcpy(fname, fprefix);
  strcat(fname, maxhfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Maximum class numbers\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{2}{c|}{$\\Delta \\equiv \\Mod{0} {4}$} &\n");
  fprintf(fp, "\\multicolumn{2}{c|}{$\\Delta \\equiv \\Mod{1} {4}$} \\\\\n");
  fprintf(fp, "\\cline{2-5}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$h$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$h$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of maximum odd parts of class numbers
  strcpy(fname, fprefix);
  strcat(fname, maxhoddfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Maximum odd parts of class numbers\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{2}{c|}{$\\Delta \\equiv \\Mod{0} {4}$} &\n");
  fprintf(fp, "\\multicolumn{2}{c|}{$\\Delta \\equiv \\Mod{1} {4}$} \\\\\n");
  fprintf(fp, "\\cline{2-5}\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$h^*$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$h^*$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of noncyclic odd part of class group
  strcpy(fname, fprefix);
  strcat(fname, ncycfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Noncyclic odd class group statistics\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$total$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$non-cyclic$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$percent$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$c(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of noncyclic odd part of class group, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, ncyc0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Noncyclic odd class group statistics, D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$total$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$non-cyclic$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$percent$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$c(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of noncyclic odd part of class group, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, ncyc1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Noncyclic odd class group statistics, D = 1(mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$total$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$non-cyclic$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$percent$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$c(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // graph points file
  strcpy(fname, fprefix);
  strcat(fname, ncycgraph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# Noncyclic odd class group statistics\n\n");
  fprintf(fp, "#x\t#c(x)\n");
  fclose(fp);

  // ********************************
  // *** 3b. p-rank probabilities ***
  // ********************************

  // file of p-rank statistics
  strcpy(fname, fprefix);
  strcat(fname, prankfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{3-9}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prank0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics, D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{3-9}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prank1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics, D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{3-9}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics (percentages)
  strcpy(fname, fprefix);
  strcat(fname, prankperfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics (percentages)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{3-9}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics (percentages), D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankper0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics (percentages), D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{3-9}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics (percentages), D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankper1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics (percentages), D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{7}{c|}{$p$} \\\\\n");
  fprintf(fp, "\\cline{3-9}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$7$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$11$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$13$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$17$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics, ratios
  strcpy(fname, fprefix);
  strcat(fname, prankratfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics (ratios)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{2,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{3,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{5,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{7,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{11,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{13,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{17,r}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics (ratios), D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankrat0file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics (ratios), D = 0 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{2,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{3,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{5,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{7,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{11,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{13,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{17,r}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of p-rank statistics (ratios), D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankrat1file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "prank statistics (ratios), D = 1 (mod 4)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r||r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$r$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{2,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{3,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{5,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{7,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{11,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{13,r}(x)$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$pr_{17,r}(x)$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank22graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 2-rank = 2\n\n");
  fprintf(fp, "#x\t#p22(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank32graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 3-rank = 2\n\n");
  fprintf(fp, "#x\t#p32(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank24graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 2-rank = 4\n\n");
  fprintf(fp, "#x\t#p24(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank52graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 5-rank = 2\n\n");
  fprintf(fp, "#x\t#p52(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank72graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 7-rank = 2\n\n");
  fprintf(fp, "#x\t#p72(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank23graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 2-rank = 3\n\n");
  fprintf(fp, "#x\t#p23(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank33graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 3-rank = 3\n\n");
  fprintf(fp, "#x\t#p33(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank53graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 5-rank = 3\n\n");
  fprintf(fp, "#x\t#p53(x)\n");
  fclose(fp);

  strcpy(fname, fprefix);
  strcat(fname, prank73graph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "# 7-rank = 3\n\n");
  fprintf(fp, "#x\t#p73(x)\n");
  fclose(fp);

  // ***************************
  // *** 4. first occurences ***
  // ***************************

  // files of first occurences of rank 2 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first22file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 2 2-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 2 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp2file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 2 p-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 3 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first23file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 3 2-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 3 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp3file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 3 p-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 4 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first24file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 4 2-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 4 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp4file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 3 p-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 5 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first25file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 5 2-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 5 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp5file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 5 p-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 6 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first26file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 6 2-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_6$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences of rank 6 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp6file);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of rank 6 p-Sylow groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_5$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$e_6$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences doubly noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, twononcycfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of doubly noncyclic groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences trebly noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, threenoncycfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of trebly noncyclic groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // files of first occurences quadruply noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, fournoncycfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "First occurences of quadruply noncyclic groups\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p_1$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_2$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_3$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$p_4$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // ****************************************
  // *** 5. number of generators required ***
  // ****************************************

  // file of maximum prime ideal required
  strcpy(fname, fprefix);
  strcat(fname, maxpfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Maximum prime ideal statistics\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\max{p}$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$local \\max{p}$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$local \\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // graph file of maximum print ideal required
  strcpy(fname, fprefix);
  strcat(fname, maxpgraph);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "#Maximum prime ideal statistics\n");
  fclose(fp);

  // file of maximum prime ideal ratios
  strcpy(fname, fprefix);
  strcat(fname, maxratfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "Maximum prime ideal statistics (ratios)\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{\\mbox{}} &\n");
  fprintf(fp, "\\multicolumn{3}{c|}{$\\max{p} / \\log{\\Delta}$} &\n");
  fprintf(fp, "\\multicolumn{3}{c|}{$\\max{p} / \\log^2 \\Delta$} \\\\\n");
  fprintf(fp, "\\cline{2-7}\n");
  fprintf(fp, "\\multicolumn{1}{|c||}{$x$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ave$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\max$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$ave$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\max$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{$\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of psplit
  strcpy(fname, fprefix);
  strcat(fname, psplitfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "First occurences of D needing maxp = p\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$p$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // file of kpsplit
  strcpy(fname, fprefix);
  strcat(fname, kpsplitfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "First occurences of D needing nump = k\n\n");
  fprintf(fp, "\\begin{table}[htb]\n");
  fprintf(fp, "\\begin{center}\n");
  fprintf(fp, "\\begin{tabular}{|r|r|r|r|r|}\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\multicolumn{1}{|c|}{$k$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# even $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{first odd $\\Delta$} &\n");
  fprintf(fp, "\\multicolumn{1}{c|}{\\# odd $\\Delta$} \\\\\n");
  fprintf(fp, "\\hline\n");
  fprintf(fp, "\\hline\n");
  fclose(fp);

  // qabc
  strcpy(fname, fprefix);
  strcat(fname, smallhfile);
  if ((fp = fopen(fname, "w")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "10 smallest h\n\n");
  fclose(fp);
}

// function to merge non-interval data files
int output_merged(int numfiles, iq_data data) {

  int i, j, key, m8;
  int e1, e2, e3, e4, e5, e6;
  long long D;
  float L1X;
  float LLI, ULI, LD, LLID, ULID;
  char fname[50], line[100];
  FILE *fp, *fp0, *fp1, *fp3, *fp4, *fp5, *fp6, *fpp, *fpp0, *fpp1, *fpp3,
      *fpp4, *fpp5, *fpp6;
  FILE *fpl0, *fpl1, *fpl5, *fpu0, *fpu1, *fpu5;
  iq_data newdata;
  float c;
  PrimeSeq sp;
  long p;

  cout << "Merging files ..." << endl;

  // ***********************
  // *** 1. L(1,X) stats ***
  // ***********************

  // successive L(1) minima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, minL0file);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive L(1) minima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minL1file);
  if ((fp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive L(1) minima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minL5file);
  if ((fp5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive L(1) maxima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, maxL0file);
  if ((fpp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive L(1) maxima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxL1file);
  if ((fpp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive L(1) maxima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxL5file);
  if ((fpp5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive LLI minima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, minLLI0file);
  if ((fpl0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive LLI minima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minLLI1file);
  if ((fpl1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive LLI minima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minLLI5file);
  if ((fpl5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive ULI maxima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, maxULI0file);
  if ((fpu0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive ULI maxima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxULI1file);
  if ((fpu1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // successive ULI maxima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxULI5file);
  if ((fpu5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // process temp files
  for (i = 0; i < numfiles; i++) {

    // open file
    sprintf(fname, "%s%d.tex", datafile, i);
    if ((fp = fopen(fname, "r")) == NULL) {
      printf("Cannot open %s.\n", fname);
      exit(1);
    }

    // process each line
    while (fgets(line, 100, fp)) {
      sscanf(line, "%d %lld %f", &key, &D, &L1X);

      // for odd D
      if (D & 1) {

        // compute LD(1)
        m8 = D % 8;
        if (m8 < 0)
          m8 += 8;
        if (m8 == 1)
          LD = L1X / 2.0;
        else
          LD = L1X * 1.5;

        // compute LLI, LLID
        c = 12 * ANTL::exp(gamma) / to<float>(ComputePi_RR() * ComputePi_RR());
        LLI = L1X * c * ANTL::log(ANTL::log(to<float>(D)));
        LLID = LD * c * ANTL::log(ANTL::log(to<float>(D << 2)));

        // compute ULI, ULID
        c = 2 * ANTL::exp(gamma);
        ULI = L1X / (c * ANTL::log(ANTL::log(to<float>(D))));
        ULID = L1X / (c * ANTL::log(ANTL::log(to<float>(D << 2))));
      }
      // for even D
      else {

        // compute LLI, ULI
        c = 8 * ANTL::exp(gamma) / to<float>(ComputePi_RR() * ComputePi_RR());
        LLI = L1X * c * ANTL::log(ANTL::log(to<float>(D)));
        c = ANTL::exp(gamma);
        ULI = L1X / (c * ANTL::log(ANTL::log(to<float>(D))));
        LLID = ULID = LD = 0;
      }

      // print out successive min / max
      switch (key) {
      case minL0_c:
        if (L1X < newdata.minL0) {
          fprintf(fp0, "%lld & %.5f & %.5f \\\\ \\hline\n", D, L1X, LLI);
          newdata.minL0 = L1X;
        }
        break;

      case minL1_c:
        if (L1X < newdata.minL1) {
          fprintf(fp1, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, LLI, LD, LLID);
          newdata.minL1 = L1X;
        }
        break;

      case minL5_c:
        if (L1X < newdata.minL5) {
          fprintf(fp5, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, LLI, LD, LLID);
          newdata.minL5 = L1X;
        }
        break;

      case maxL0_c:
        if (L1X > newdata.maxL0) {
          fprintf(fpp0, "%lld & %.5f & %.5f \\\\ \\hline\n", D, L1X, ULI);
          newdata.maxL0 = L1X;
        }
        break;

      case maxL1_c:
        if (L1X > newdata.maxL1) {
          fprintf(fpp1, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, ULI, LD, ULID);
          newdata.maxL1 = L1X;
        }
        break;

      case maxL5_c:
        if (L1X > newdata.maxL5) {
          fprintf(fpp5, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, ULI, LD, ULID);
          newdata.maxL5 = L1X;
        }
        break;

      case minLLI0_c:
        if (LLI < newdata.minLLI0) {
          fprintf(fpl0, "%lld & %.5f & %.5f \\\\ \\hline\n", D, L1X, LLI);
          newdata.minLLI0 = LLI;
        }
        break;

      case minLLI1_c:
        if (LLI < newdata.minLLI1) {
          fprintf(fpl1, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, LLI, LD, LLID);
          newdata.minLLI1 = LLI;
        }
        break;

      case minLLI5_c:
        if (LLI < newdata.minLLI5) {
          fprintf(fpl5, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, LLI, LD, LLID);
          newdata.minLLI5 = LLI;
        }
        break;

      case maxULI0_c:
        if (ULI > newdata.maxULI0) {
          fprintf(fpu0, "%lld & %.5f & %.5f \\\\ \\hline\n", D, L1X, ULI);
          newdata.maxULI0 = ULI;
        }
        break;

      case maxULI1_c:
        if (ULI > newdata.maxULI1) {
          fprintf(fpu1, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, ULI, LD, ULID);
          newdata.maxULI1 = ULI;
        }
        break;

      case maxULI5_c:
        if (ULI > newdata.maxULI5) {
          fprintf(fpu5, "%lld & %.5f & %.5f & %.5f & %.5f \\\\ \\hline\n", D,
                  L1X, ULI, LD, ULID);
          newdata.maxULI5 = ULI;
        }
        break;

      default:
        cerr << "Error reading file -- invalid key" << endl;
        break;
      }
    }
    fclose(fp);
  }

  fclose(fp0);
  fclose(fp1);
  fclose(fp5);
  fclose(fpp0);
  fclose(fpp1);
  fclose(fpp5);
  fclose(fpl0);
  fclose(fpl1);
  fclose(fpl5);
  fclose(fpu0);
  fclose(fpu1);
  fclose(fpu5);

  // ***************************
  // *** 4. First occurences ***
  // ***************************

  // files of first occurences of rank 2 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first22file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 2 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp2file);
  if ((fpp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 3 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first23file);
  if ((fp3 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 3 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp3file);
  if ((fpp3 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 4 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first24file);
  if ((fp4 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 4 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp4file);
  if ((fpp4 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 5 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first25file);
  if ((fp5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 5 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp5file);
  if ((fpp5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 6 2-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, first26file);
  if ((fp6 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences of rank 6 p-Sylow subgroups
  strcpy(fname, fprefix);
  strcat(fname, firstp6file);
  if ((fpp6 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print out 2-Sylow group statistics
  for (e1 = 1; e1 < MAXE1; e1++)
    for (e2 = 1; e2 < MAXE2; e2++) {
      // rank 2 2-Sylow subgroups
      if (data.prank0[0][e1][e2][0][0][0][0] > 0 ||
          data.prank1[0][e1][e2][0][0][0][0] > 0) {
        fprintf(fp, "%d & %d", e1, e2);

        if (data.prank0[0][e1][e2][0][0][0][0] > 0) {
          fprintf(fp, " & %lld & %lld",
                  to<long long>(data.first_prank0[0][e1][e2][0][0][0][0]),
                  data.prank0[0][e1][e2][0][0][0][0]);
        } else {
          fprintf(fp, " & * & *");
        }
        if (data.prank1[0][e1][e2][0][0][0][0] > 0) {
          fprintf(fp, " & %lld & %lld \\\\ \\hline\n",
                  to<long long>(data.first_prank1[0][e1][e2][0][0][0][0]),
                  data.prank1[0][e1][e2][0][0][0][0]);
        } else {
          fprintf(fp, " & * & * \\\\ \\hline\n");
        }
      }

      for (e3 = 1; e3 < MAXE3; e3++) {
        // rank 3 2-Sylow subgroups
        if (data.prank0[0][e1][e2][e3][0][0][0] > 0 ||
            data.prank1[0][e1][e2][e3][0][0][0] > 0) {
          fprintf(fp3, "%d & %d & %d", e1, e2, e3);

          if (data.prank0[0][e1][e2][e3][0][0][0] > 0) {
            fprintf(fp3, " & %lld & %lld",
                    to<long long>(data.first_prank0[0][e1][e2][e3][0][0][0]),
                    data.prank0[0][e1][e2][e3][0][0][0]);
          } else {
            fprintf(fp3, " & * & *");
          }
          if (data.prank1[0][e1][e2][e3][0][0][0] > 0) {
            fprintf(fp3, " & %lld & %lld \\\\ \\hline\n",
                    to<long long>(data.first_prank1[0][e1][e2][e3][0][0][0]),
                    data.prank1[0][e1][e2][e3][0][0][0]);
          } else {
            fprintf(fp3, " & * & * \\\\ \\hline\n");
          }
        }
        for (e4 = 1; e4 < MAXE4; e4++) {
          // rank 4 2-Sylow groups
          if (data.prank0[0][e1][e2][e3][e4][0][0] > 0 ||
              data.prank1[0][e1][e2][e3][e4][0][0] > 0) {
            fprintf(fp4, "%d & %d & %d & %d", e1, e2, e3, e4);

            if (data.prank0[0][e1][e2][e3][e4][0][0] > 0) {
              fprintf(fp4, " & %lld & %lld",
                      to<long long>(data.first_prank0[0][e1][e2][e3][e4][0][0]),
                      data.prank0[0][e1][e2][e3][e4][0][0]);
            } else {
              fprintf(fp4, " & * & *");
            }
            if (data.prank1[0][e1][e2][e3][e4][0][0] > 0) {
              fprintf(fp4, " & %lld & %lld \\\\ \\hline\n",
                      to<long long>(data.first_prank1[0][e1][e2][e3][e4][0][0]),
                      data.prank1[0][e1][e2][e3][e4][0][0]);
            } else {
              fprintf(fp4, " & * & * \\\\ \\hline\n");
            }
          }
          for (e5 = 1; e5 < MAXE5; e5++) {
            // rank 5 2-Sylow groups
            if (data.prank0[0][e1][e2][e3][e4][e5][0] > 0 ||
                data.prank1[0][e1][e2][e3][e4][e5][0] > 0) {
              fprintf(fp5, "%d & %d & %d & %d & %d", e1, e2, e3, e4, e5);

              if (data.prank0[0][e1][e2][e3][e4][e5][0] > 0) {
                fprintf(
                    fp5, " & %lld & %lld",
                    to<long long>(data.first_prank0[0][e1][e2][e3][e4][e5][0]),
                    data.prank0[0][e1][e2][e3][e4][e5][0]);
              } else {
                fprintf(fp5, " & * & *");
              }
              if (data.prank1[0][e1][e2][e3][e4][e5][0] > 0) {
                fprintf(
                    fp5, " & %lld & %lld \\\\ \\hline\n",
                    to<long long>(data.first_prank1[0][e1][e2][e3][e4][e5][0]),
                    data.prank1[0][e1][e2][e3][e4][e5][0]);
              } else {
                fprintf(fp5, " & * & * \\\\ \\hline\n");
              }
            }
            for (e6 = 1; e6 < MAXE6; e6++) {
              // rank 6 2-Sylow groups
              if (data.prank0[0][e1][e2][e3][e4][e5][e6] > 0 ||
                  data.prank1[0][e1][e2][e3][e4][e5][e6] > 0) {
                fprintf(fp6, "%d & %d & %d & %d & %d & %d", e1, e2, e3, e4, e5,
                        e6);

                if (data.prank0[0][e1][e2][e3][e4][e5][e6] > 0) {
                  fprintf(fp6, " & %lld & %lld",
                          to<long long>(
                              data.first_prank0[0][e1][e2][e3][e4][e5][e6]),
                          data.prank0[0][e1][e2][e3][e4][e5][e6]);
                } else {
                  fprintf(fp6, " & * & *");
                }
                if (data.prank1[0][e1][e2][e3][e4][e5][e6] > 0) {
                  fprintf(fp6, " & %lld & %lld \\\\ \\hline\n",
                          to<long long>(
                              data.first_prank1[0][e1][e2][e3][e4][e5][e6]),
                          data.prank1[0][e1][e2][e3][e4][e5][e6]);
                } else {
                  fprintf(fp6, " & * & * \\\\ \\hline\n");
                }
              }
            }
          }
        }
      }
    }

  // print out p-Sylow group statistics
  sp.reset(3);
  for (j = 1; j < MAXPLARGE; j++) {
    p = sp.next();
    for (e1 = 1; e1 < MAXE1; e1++)
      for (e2 = 1; e2 < MAXE2; e2++) {
        // p - Sylow groups
        if (data.prank0[j][e1][e2][0][0][0][0] > 0 ||
            data.prank1[j][e1][e2][0][0][0][0] > 0) {
          fprintf(fpp, "%d & %d & %d", p, e1, e2);

          if (data.prank0[j][e1][e2][0][0][0][0] > 0) {
            fprintf(fpp, " & %lld & %lld",
                    to<long long>(data.first_prank0[j][e1][e2][0][0][0][0]),
                    data.prank0[j][e1][e2][0][0][0][0]);
          } else {
            fprintf(fpp, " & * & *");
          }
          if (data.prank1[j][e1][e2][0][0][0][0] > 0) {
            fprintf(fpp, " & %lld & %lld \\\\ \\hline\n",
                    to<long long>(data.first_prank1[j][e1][e2][0][0][0][0]),
                    data.prank1[j][e1][e2][0][0][0][0]);
          } else {
            fprintf(fpp, " & * & * \\\\ \\hline\n");
          }
        }
        for (e3 = 1; e3 < MAXE3; e3++) {
          // p-Sylow groups
          if (data.prank0[j][e1][e2][e3][0][0][0] > 0 ||
              data.prank1[j][e1][e2][e3][0][0][0] > 0) {
            fprintf(fpp3, "%d & %d & %d & %d", p, e1, e2, e3);

            if (data.prank0[j][e1][e2][e3][0][0][0] > 0) {
              fprintf(fpp3, " & %lld & %lld",
                      to<long long>(data.first_prank0[j][e1][e2][e3][0][0][0]),
                      data.prank0[j][e1][e2][e3][0][0][0]);
            } else {
              fprintf(fpp3, " & * & *");
            }
            if (data.prank1[j][e1][e2][e3][0][0][0] > 0) {
              fprintf(fpp3, " & %lld & %lld \\\\ \\hline\n",
                      to<long long>(data.first_prank1[j][e1][e2][e3][0][0][0]),
                      data.prank1[j][e1][e2][e3][0][0][0]);
            } else {
              fprintf(fpp3, " & * & * \\\\ \\hline\n");
            }
          }
          for (e4 = 1; e4 < MAXE4; e4++) {
            // p-Sylow groups
            if (data.prank0[j][e1][e2][e3][e4][0][0] > 0 ||
                data.prank1[j][e1][e2][e3][e4][0][0] > 0) {
              fprintf(fpp4, "%d & %d & %d & %d & %d", p, e1, e2, e3, e4);

              if (data.prank0[j][e1][e2][e3][e4][0][0] > 0) {
                fprintf(
                    fpp4, " & %lld & %lld",
                    to<long long>(data.first_prank0[j][e1][e2][e3][e4][0][0]),
                    data.prank0[j][e1][e2][e3][e4][0][0]);
              } else {
                fprintf(fpp4, " & * & *");
              }
              if (data.prank1[j][e1][e2][e3][e4][0][0] > 0) {
                fprintf(
                    fpp4, " & %lld & %lld \\\\ \\hline\n",
                    to<long long>(data.first_prank1[j][e1][e2][e3][e4][0][0]),
                    data.prank1[j][e1][e2][e3][e4][0][0]);
              } else {
                fprintf(fpp4, " & * & * \\\\ \\hline\n");
              }
            }
            for (e5 = 1; e5 < MAXE5; e5++) {
              // p-Sylow groups
              if (data.prank0[j][e1][e2][e3][e4][e5][0] > 0 ||
                  data.prank1[j][e1][e2][e3][e4][e5][0] > 0) {
                fprintf(fpp5, "%d & %d & %d & %d & %d & %d", p, e1, e2, e3, e4,
                        e5);

                if (data.prank0[j][e1][e2][e3][e4][e5][0] > 0) {
                  fprintf(fpp5, " & %lld & %lld",
                          to<long long>(
                              data.first_prank0[j][e1][e2][e3][e4][e5][0]),
                          data.prank0[j][e1][e2][e3][e4][e5][0]);
                } else {
                  fprintf(fpp5, " & * & *");
                }
                if (data.prank1[j][e1][e2][e3][e4][e5][0] > 0) {
                  fprintf(fpp5, " & %lld & %lld \\\\ \\hline\n",
                          to<long long>(
                              data.first_prank1[j][e1][e2][e3][e4][e5][0]),
                          data.prank1[j][e1][e2][e3][e4][e5][0]);
                } else {
                  fprintf(fpp5, " & * & * \\\\ \\hline\n");
                }
              }
              for (e6 = 1; e6 < MAXE6; e6++) {
                // p-Sylow groups
                if (data.prank0[j][e1][e2][e3][e4][e5][e6] > 0 ||
                    data.prank1[j][e1][e2][e3][e4][e5][e6] > 0) {
                  fprintf(fpp6, "%d & %d & %d & %d & %d & %d & %d", p, e1, e2,
                          e3, e4, e5, e6);

                  if (data.prank0[j][e1][e2][e3][e4][e5][e6] > 0) {
                    fprintf(fpp6, " & %lld & %lld",
                            to<long long>(
                                data.first_prank0[j][e1][e2][e3][e4][e5][e6]),
                            data.prank0[j][e1][e2][e3][e4][e5][e6]);
                  } else {
                    fprintf(fpp6, " & * & *");
                  }
                  if (data.prank1[j][e1][e2][e3][e4][e5][e6] > 0) {
                    fprintf(fpp6, " & %lld & %lld \\\\ \\hline\n",
                            to<long long>(
                                data.first_prank1[j][e1][e2][e3][e4][e5][e6]),
                            data.prank1[j][e1][e2][e3][e4][e5][e6]);
                  } else {
                    fprintf(fpp6, " & * & * \\\\ \\hline\n");
                  }
                }
              }
            }
          }
        }
      }
  }
  fclose(fp);
  fclose(fpp);
  fclose(fp3);
  fclose(fpp3);
  fclose(fp4);
  fclose(fpp4);
  fclose(fp5);
  fclose(fpp5);
  fclose(fp6);
  fclose(fpp6);

  // files of first occurences doubly noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, twononcycfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences trebly noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, threenoncycfile);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // files of first occurences quadruply noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, fournoncycfile);
  if ((fp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  long primeslist[LASTPRIME];
  sp.reset(2);
  for (j = 0; j < LASTPRIME; j++) {
    p = sp.next();
    primeslist[j] = p;
  }

  for (e1 = 0; e1 < MAXP1; e1++)
    for (e2 = 0; e2 < MAXP2; e2++)
      for (e3 = 0; e3 < MAXP3; e3++)
        for (e4 = 0; e4 < MAXP4; e4++) {

          if (e3 == 0) {
            // doubly non-cyclic
            if (data.two_noncyc0[e1][e2][e3][e4] > 0 ||
                data.two_noncyc1[e1][e2][e3][e4] > 0) {
              fprintf(fp, "%ld & %ld & ", primeslist[e1], primeslist[e2]);

              if (data.two_noncyc0[e1][e2][e3][e4] > 0) {
                fprintf(fp, "%lld & %lld & ",
                        to<long long>(data.first_two_noncyc0[e1][e2][e3][e4]),
                        data.two_noncyc0[e1][e2][e3][e4]);
              } else {
                fprintf(fp, " * & * & ");
              }

              if (data.two_noncyc1[e1][e2][e3][e4] > 0) {
                fprintf(fp, "%lld & %lld & \\\\ \\hline\n",
                        to<long long>(data.first_two_noncyc1[e1][e2][e3][e4]),
                        data.two_noncyc1[e1][e2][e3][e4]);
              } else {
                fprintf(fp, " * & * & \\\\ \\hline\n");
              }
            }
          } else if (e4 == 0) {
            // trebly non-cyclic
            if (data.two_noncyc0[e1][e2][e3][e4] > 0 ||
                data.two_noncyc1[e1][e2][e3][e4] > 0) {
              fprintf(fp0, "%ld & %ld & %ld & ", primeslist[e1], primeslist[e2],
                      primeslist[e3]);

              if (data.two_noncyc0[e1][e2][e3][e4] > 0) {
                fprintf(fp0, "%lld & %lld & ",
                        to<long long>(data.first_two_noncyc0[e1][e2][e3][e4]),
                        data.two_noncyc0[e1][e2][e3][e4]);
              } else {
                fprintf(fp0, " * & * & ");
              }

              if (data.two_noncyc1[e1][e2][e3][e4] > 0) {
                fprintf(fp0, "%lld & %lld & \\\\ \\hline\n",
                        to<long long>(data.first_two_noncyc1[e1][e2][e3][e4]),
                        data.two_noncyc1[e1][e2][e3][e4]);
              } else {
                fprintf(fp0, " * & * & \\\\ \\hline\n");
              }
            }
          } else {
            // quadruply non-cyclic
            if (data.two_noncyc0[e1][e2][e3][e4] > 0 ||
                data.two_noncyc1[e1][e2][e3][e4] > 0) {
              fprintf(fp1, "%ld & %ld & %ld & %ld & ", primeslist[e1],
                      primeslist[e2], primeslist[e3], primeslist[e4]);

              if (data.two_noncyc0[e1][e2][e3][e4] > 0) {
                fprintf(fp1, "%lld & %lld & ",
                        to<long long>(data.first_two_noncyc0[e1][e2][e3][e4]),
                        data.two_noncyc0[e1][e2][e3][e4]);
              } else {
                fprintf(fp1, " * & * & ");
              }

              if (data.two_noncyc1[e1][e2][e3][e4] > 0) {
                fprintf(fp1, "%lld & %lld & \\\\ \\hline\n",
                        to<long long>(data.first_two_noncyc1[e1][e2][e3][e4]),
                        data.two_noncyc1[e1][e2][e3][e4]);
              } else {
                fprintf(fp1, " * & * & \\\\ \\hline\n");
              }
            }
          }
        }

  fclose(fp);
  fclose(fp0);
  fclose(fp1);

  // 5. Generators

  // files of first occurences of p splits
  strcpy(fname, fprefix);
  strcat(fname, psplitfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // files of first occurences of k p splits
  strcpy(fname, fprefix);
  strcat(fname, kpsplitfile);
  if ((fpp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  sp.reset(3);
  for (j = 0; j < MAXPLARGE; j++) {
    p = sp.next();
    fprintf(fp, "%d & ", p);
    if (data.psplit0[j] > 0)
      fprintf(fp, "%lld & %lld & ", to<long long>(data.first_psplit0[j]),
              data.psplit0[j]);
    else
      fprintf(fp, "* & * & ");
    if (data.psplit1[j] > 0)
      fprintf(fp, "%lld & %lld \\\\ \\hline\n",
              to<long long>(data.first_psplit1[j]), data.psplit1[j]);
    else
      fprintf(fp, "* & * \\\\ \\hline\n");

    fprintf(fpp, "%d & ", j);
    if (data.kpsplit0[j] > 0)
      fprintf(fpp, "%lld & %lld & ", to<long long>(data.first_kpsplit0[j]),
              data.kpsplit0[j]);
    else
      fprintf(fpp, "* & * & ");
    if (data.kpsplit1[j] > 0)
      fprintf(fpp, "%lld & %lld \\\\ \\hline\n",
              to<long long>(data.first_kpsplit1[j]), data.kpsplit1[j]);
    else
      fprintf(fpp, "* & * \\\\ \\hline\n");
  }

  return 0;
}

// function to output interval data to latex files
int output_interval(iq_data data, ZZ thresh) {

  // parameters
  char fname[50];
  FILE *fp, *fp0, *fp1, *fp5, *fpp, *fpp0, *fpp1, *fpr, *fpr0, *fpr1;
  FILE *fps, *fps0, *fps1, *fpc, *fpc0, *fpc1;
  double per, per0, per1;
  double rat, rat0, rat1;
  float maxULI, minLLI;
  ZZ minD, maxD;
  PrimeSeq sp;
  long p;
  double prob[MAXP][3];
  int i, j;

  // **************************
  // *** 1. L(1,X) function ***
  // **************************

  // file of average L function
  strcpy(fname, fprefix);
  strcat(fname, aveLfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print out statistics
  fprintf(fp, "%lld", to<long long>(thresh));
  per =
      to<double>((data.sumL0 + data.sumL1) / to<RR>(data.total0 + data.total1));
  per0 = to<double>((data.sumL0) / to<RR>(data.total0));
  per1 = to<double>((data.sumL1) / to<RR>(data.total1));
  fprintf(fp, " & %.5f & %.5f & %.5f \\\\\n", per, per0, per1);
  fclose(fp);

  // ***********************************
  // *** 2. divisors of h statistics ***
  // ***********************************

  // divisors of h
  strcpy(fname, fprefix);
  strcat(fname, pdivsfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // divisors of h, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivs0file);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivs1file);
  if ((fp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h percentages
  strcpy(fname, fprefix);
  strcat(fname, pdivsperfile);
  if ((fpp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h percentages, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivsper0file);
  if ((fpp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h percentages, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivsper1file);
  if ((fpp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivsratfile);
  if ((fpr = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivsrat0file);
  if ((fpr0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivsrat1file);
  if ((fpr1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, (squares)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqrfile);
  if ((fps = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, (squares) D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqr0file);
  if ((fps0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, (squares) D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqr1file);
  if ((fps1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, (cubes)
  strcpy(fname, fprefix);
  strcat(fname, pdivscubfile);
  if ((fpc = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, (cubes) D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivscub0file);
  if ((fpc0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // divisors of h ratios, (cubes) D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivscub1file);
  if ((fpc1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s. \n", fname);
    exit(1);
  }
  // calculate pr(p|h)
  sp.reset(3);
  p = sp.next();
  for (i = 1; i < MAXP; i++) {
    prob[i][0] = 100.0 * (1 - nu(p));
    p = sp.next();
  }

  // calculate pr(p^2|h)
  sp.reset(3);
  p = sp.next();
  for (i = 1; i < MAXP; i++) {
    prob[i][1] = p * nu(p);
    prob[i][1] /= (p - 1);
    prob[i][1] = 1 - prob[i][1];
    p = sp.next();
  }

  // calculate pr(p^3|h)
  sp.reset(3);
  p = sp.next();
  for (i = 1; i < MAXP; i++) {
    prob[i][2] = p * p * p * nu(p);
    prob[i][2] /= ((p - 1) * (p - 1) * (p + 1));
    prob[i][2] = 1 - prob[i][2];
    p = sp.next();
  }

  // print out statistics
  fprintf(fp, "%lld", to<long long>(thresh));
  fprintf(fp0, "%lld", to<long long>(thresh));
  fprintf(fp1, "%lld", to<long long>(thresh));
  fprintf(fpp, "%lld", to<long long>(thresh));
  fprintf(fpp0, "%lld", to<long long>(thresh));
  fprintf(fpp1, "%lld", to<long long>(thresh));
  fprintf(fpr, "%lld", to<long long>(thresh));
  fprintf(fpr0, "%lld", to<long long>(thresh));
  fprintf(fpr1, "%lld", to<long long>(thresh));
  fprintf(fps, "%lld", to<long long>(thresh));
  fprintf(fps0, "%lld", to<long long>(thresh));
  fprintf(fps1, "%lld", to<long long>(thresh));
  fprintf(fpc, "%lld", to<long long>(thresh));
  fprintf(fpc0, "%lld", to<long long>(thresh));
  fprintf(fpc1, "%lld", to<long long>(thresh));
  for (i = 1; i < MAXP1; i++) {
    // divisors of h
    fprintf(fp, " & %lld", data.pdivs0[i][0] + data.pdivs1[i][0]);
    fprintf(fp0, " & %lld", data.pdivs0[i][0]);
    fprintf(fp1, " & %lld", data.pdivs1[i][0]);
    // divisors of h percentages
    per = 100.0 * to<double>(to<RR>(data.pdivs0[i][0] + data.pdivs1[i][0]) /
                             to<RR>(data.total0 + data.total1));
    per0 = 100.0 * to<double>(to<RR>(data.pdivs0[i][0]) / to<RR>(data.total0));
    per1 = 100.0 * to<double>(to<RR>(data.pdivs1[i][0]) / to<RR>(data.total1));
    fprintf(fpp, " & %.5f", per);
    fprintf(fpp0, " & %.5f", per0);
    fprintf(fpp1, " & %.5f", per1);
    // divisors of h ratios
    rat = to<double>(per / prob[i][0]);
    rat0 = to<double>(per0 / prob[i][0]);
    rat1 = to<double>(per1 / prob[i][0]);
    fprintf(fpr, " & %.5f", rat);
    fprintf(fpr0, " & %.5f", rat0);
    fprintf(fpr1, " & %.5f", rat1);
    // divisors of h ratios (squares)
    fprintf(fps, " & %.5f",
            to<double>(to<RR>(data.pdivs0[i][1] + data.pdivs1[i][1]) /
                       (prob[i][1] * to<RR>(data.total0 + data.total1))));
    fprintf(fps0, " & %.5f",
            to<double>(to<RR>(data.pdivs0[i][1]) /
                       (prob[i][1] * to<RR>(data.total0))));
    fprintf(fps1, " & %.5f",
            to<double>(to<RR>(data.pdivs1[i][1]) /
                       (prob[i][1] * to<RR>(data.total1))));
    // divisors of h ratios (cubes)
    fprintf(fpc, " & %.5f",
            to<double>(to<RR>(data.pdivs0[i][2] + data.pdivs1[i][2]) /
                       (prob[i][2] * to<RR>(data.total0 + data.total1))));
    fprintf(fpc0, " & %.5f",
            to<double>(to<RR>(data.pdivs0[i][2]) /
                       (prob[i][2] * to<RR>(data.total0))));
    fprintf(fpc1, " & %.5f",
            to<double>(to<RR>(data.pdivs1[i][2]) /
                       (prob[i][2] * to<RR>(data.total1))));
  }
  // close files
  fprintf(fp, " \\\\ \\hline\n");
  fprintf(fp0, " \\\\ \\hline\n");
  fprintf(fp1, " \\\\ \\hline\n");
  fprintf(fpp, " \\\\ \\hline\n");
  fprintf(fpp0, " \\\\ \\hline\n");
  fprintf(fpp1, " \\\\ \\hline\n");
  fprintf(fpr, " \\\\ \\hline\n");
  fprintf(fpr0, " \\\\ \\hline\n");
  fprintf(fpr1, " \\\\ \\hline\n");
  fprintf(fps, "\\\\ \\hline\n");
  fprintf(fps0, " \\\\ \\hline\n");
  fprintf(fps1, " \\\\ \\hline\n");
  fprintf(fpc, " \\\\ \\hline\n");
  fprintf(fpc0, " \\\\ \\hline\n");
  fprintf(fpc1, " \\\\ \\hline\n");
  fclose(fp);
  fclose(fp0);
  fclose(fp1);
  fclose(fpp);
  fclose(fpp0);
  fclose(fpp1);
  fclose(fpr);
  fclose(fpr0);
  fclose(fpr1);
  fclose(fps);
  fclose(fps0);
  fclose(fps1);
  fclose(fpc);
  fclose(fpc0);
  fclose(fpc1);

  // **********************************************
  // *** 3a. odd part of class group statistics ***
  // **********************************************

  // file of maximum class numbers
  strcpy(fname, fprefix);
  strcat(fname, maxhfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of maximum odd parts of class numbers
  strcpy(fname, fprefix);
  strcat(fname, maxhoddfile);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print our statistics
  fprintf(fp, "%lld & %ld & %lld & %ld & %lld \\\\ \\hline\n",
          to<long long>(thresh), to<long>(data.maxh0),
          to<long long>(data.maxh0D), to<long>(data.maxh1),
          to<long long>(data.maxh1D));
  fprintf(fp0, "%lld & %ld & %lld & %ld & %lld \\\\ \\hline\n",
          to<long long>(thresh), to<long>(data.maxhodd0),
          to<long long>(data.maxhodd0D), to<long>(data.maxhodd1),
          to<long>(data.maxhodd1D));

  // close files
  fclose(fp);
  fclose(fp0);

  // file of noncyclic odd part of class group
  strcpy(fname, fprefix);
  strcat(fname, ncycfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of noncyclic odd part of class group
  strcpy(fname, fprefix);
  strcat(fname, ncyc0file);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of noncyclic odd part of class group
  strcpy(fname, fprefix);
  strcat(fname, ncyc1file);
  if ((fp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print out statistics
  fprintf(fp, "%lld", to<long long>(thresh));
  fprintf(fp0, "%lld", to<long long>(thresh));
  fprintf(fp1, "%lld", to<long long>(thresh));
  // total, non-cyclic, c(x)
  per = 100.0 * to<double>(to<RR>(data.noncyc0 + data.noncyc1) /
                           to<RR>(data.total0 + data.total1));
  per0 = 100.0 * to<double>(to<RR>(data.noncyc0) / to<RR>(data.total0));
  per1 = 100.0 * to<double>(to<RR>(data.noncyc1) / to<RR>(data.total1));
  fprintf(fp, " & %lld & %lld & %.5f & %.5f \\\\ \\hline\n",
          data.total0 + data.total1, data.noncyc0 + data.noncyc1, per,
          (100 - per) / (100.0 * cycprob));
  fprintf(fp0, " & %lld & %lld & %.5f & %.5f \\\\ \\hline\n", data.total0,
          data.noncyc0, per0, (100 - per0) / (100.0 * cycprob));
  fprintf(fp1, " & %lld & %lld & %.5f & %.5f \\\\ \\hline\n", data.total1,
          data.noncyc1, per1, (100 - per1) / (100.0 * cycprob));
  // close files
  fclose(fp);
  fclose(fp0);
  fclose(fp1);

  // *****************************
  // *** 3b. p-rank statistics ***
  // *****************************

  // file of p-rank statistics
  strcpy(fname, fprefix);
  strcat(fname, prankfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics
  strcpy(fname, fprefix);
  strcat(fname, prank0file);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics
  strcpy(fname, fprefix);
  strcat(fname, prank1file);
  if ((fp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics (percentages)
  strcpy(fname, fprefix);
  strcat(fname, prankperfile);
  if ((fpp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics (percentages), D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankper0file);
  if ((fpp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics (percentages), D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankper1file);
  if ((fpp1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics (ratios)
  strcpy(fname, fprefix);
  strcat(fname, prankratfile);
  if ((fpr = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics (ratios), D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankrat0file);
  if ((fpr0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // file of p-rank statistics (ratios), D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankrat1file);
  if ((fpr1 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print out statistics
  fprintf(fp, "%lld", to<long long>(thresh));
  fprintf(fp0, "%lld", to<long long>(thresh));
  fprintf(fp1, "%lld", to<long long>(thresh));
  fprintf(fpp, "%lld", to<long long>(thresh));
  fprintf(fpp0, "%lld", to<long long>(thresh));
  fprintf(fpp1, "%lld", to<long long>(thresh));
  fprintf(fpr, "%lld", to<long long>(thresh));
  fprintf(fpr0, "%lld", to<long long>(thresh));
  fprintf(fpr1, "%lld", to<long long>(thresh));

  // check if p-rank = 2, 3, ... , 6
  long long pr0, pr1;
  int e1, e2, e3, e4, e5, e6;
  int s3, s4, s5, s6;
  int f3, f4, f5, f6;
  double ratio[MAXP], value;

  for (i = 2; i <= 6; i++) {
    fprintf(fp, " & %d", i);
    fprintf(fp0, " & %d", i);
    fprintf(fp1, " & %d", i);
    fprintf(fpp, " & %d", i);
    fprintf(fpp0, " & %d", i);
    fprintf(fpp1, " & %d", i);
    fprintf(fpr, " & %d", i);
    fprintf(fpr0, " & %d", i);
    fprintf(fpr1, " & %d", i);

    sp.reset(2);
    for (j = 0; j < MAXPSMALL; j++) {
      p = sp.next();

      // calculate ratio
      ratio[j] = nu(p);
      ratio[j] /= pow(to<double>(p), i * i);
      value = nu(p, i);
      ratio[j] /= (value * value);

      // find exponents
      pr0 = pr1 = 0;
      if (i == 2) {
        s3 = s4 = s5 = s6 = 0;
        f3 = f4 = f5 = f6 = 1;

      } else if (i == 3) {
        f3 = MAXE3;
        s3 = f4 = f5 = f6 = 1;
        s4 = s5 = s6 = 0;
      } else if (i == 4) {
        f3 = MAXE3;
        f4 = MAXE4;
        s3 = s4 = f5 = f6 = 1;
        s5 = s6 = 0;
      } else if (i == 5) {
        f3 = MAXE3;
        f4 = MAXE4;
        f5 = MAXE5;
        s3 = s4 = s5 = f6 = 1;
        s6 = 0;
      } else {
        f3 = MAXE3;
        f4 = MAXE4;
        f5 = MAXE5;
        f6 = MAXE6;
        s3 = s4 = s5 = s6 = 1;
      }

      for (e1 = 1; e1 < MAXE1; e1++)
        for (e2 = 1; e2 < MAXE2; e2++)
          for (e3 = s3; e3 < f3; e3++)
            for (e4 = s4; e4 < f4; e4++)
              for (e5 = s5; e5 < f5; e5++)
                for (e6 = s6; e6 < f6; e6++) {
                  pr0 += data.prank0[j][e1][e2][e3][e4][e5][e6];
                  pr1 += data.prank1[j][e1][e2][e3][e4][e5][e6];
                }
      // p-rank values
      fprintf(fp, " & %lld", pr0 + pr1);
      fprintf(fp0, " & %lld", pr0);
      fprintf(fp1, " & %lld", pr1);
      // percentages
      per = 100.0 *
            to<double>(to<RR>(pr0 + pr1) / to<RR>(data.total0 + data.total1));
      per0 = 100.0 * to<double>(to<RR>(pr0) / to<RR>(data.total0));
      per1 = 100.0 * to<double>(to<RR>(pr1) / to<RR>(data.total1));
      fprintf(fpp, " & %.5f", per);
      fprintf(fpp0, " & %.5f", per0);
      fprintf(fpp1, " & %.5f", per1);
      // ratios
      fprintf(fpr, " & %.5f", per / (100.0 * ratio[j]));
      fprintf(fpr0, " & %.5f", per / (100.0 * ratio[j]));
      fprintf(fpr1, " & %.5f", per / (100.0 * ratio[j]));
    }
    fprintf(fp, "\\\\ \\hline\n");
    fprintf(fp0, " \\\\ \\hline\n");
    fprintf(fp1, " \\\\ \\hline\n");
    fprintf(fpp, " \\\\ \\hline\n");
    fprintf(fpp0, " \\\\ \\hline\n");
    fprintf(fpp1, " \\\\ \\hline\n");
    fprintf(fpr, " \\\\ \\hline\n");
    fprintf(fpr0, " \\\\ \\hline\n");
    fprintf(fpr1, " \\\\ \\hline\n");
  }

  // close files
  fclose(fp);
  fclose(fp0);
  fclose(fp1);
  fclose(fpp);
  fclose(fpp0);
  fclose(fpp1);
  fclose(fpr);
  fclose(fpr0);
  fclose(fpr1);

  // *******************************
  // *** 5. number of generators ***
  // *******************************

  // file of maximum prime ideal required
  strcpy(fname, fprefix);
  strcat(fname, maxpfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // max log (and log log) ratios file
  strcpy(fname, fprefix);
  strcat(fname, maxratfile);
  if ((fp0 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // print out statistics
  fprintf(fp, "%lld & %ld & %lld & %ld & %lld \\\\ \\hline\n",
          to<long long>(thresh), to<long>(data.maxp), to<long long>(data.maxpD),
          to<long>(data.localmaxp), to<long long>(data.localmaxpD));
  fprintf(
      fp0, "%lld & %.5f & %.5f & %lld", to<long long>(thresh),
      to<double>(to<RR>(data.sumlograt) / to<RR>(data.total0 + data.total1)),
      to<double>(data.maxlograt), to<long long>(data.maxlogratD));
  fprintf(
      fp0, " & %.5f & %.5f & %lld \\\\ \\hline\n",
      to<double>(to<RR>(data.sumloglograt) / to<RR>(data.total0 + data.total1)),
      to<double>(data.maxloglograt), to<long long>(data.maxloglogratD));
  fclose(fp);
  fclose(fp0);

  // qabc
  strcpy(fname, fprefix);
  strcat(fname, smallhfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "%lld & ", to<long long>(thresh));
  for (i = 0; i < NUMH; i++) {
    fprintf(fp, "%ld & %lld & ", data.smallh[i],
            to<long long>(data.smallhD[i]));
  }
  fprintf(fp, "\n");
  fclose(fp);

  return 0;
}

// function to output interval data to graph files
int output_graph_interval(iq_data data, ZZ gthresh) {

  char fname[50];
  FILE *fp, *fp3, *fp5, *fp7;
  PrimeSeq sp;
  long p;
  double prob[MAXP];
  int i, j;
  int e1, e2, e3, e4, e5, e6;
  int s3, s4, s5, s6, f3, f4, f5, f6;
  double ratio[MAXP], value;
  long long pr0, pr1;
  float maxULI, minLLI;
  ZZ minD, maxD;

  // file of min LLI
  strcpy(fname, fprefix);
  strcat(fname, localLLIgraph);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print out statistics
  minLLI = data.minLLI0;
  minD = data.minLLI0D;
  if (data.minLLI1 < minLLI) {
    minLLI = data.minLLI1;
    minD = data.minLLI1D;
  }
  if (data.minLLI5 < minLLI) {
    minLLI = data.minLLI5;
    minD = data.minLLI5D;
  }
  fprintf(fp, "%lld & %.5f \\\\ \n", to<long long>(minD), minLLI);
  fclose(fp);

  // file of max ULI
  strcpy(fname, fprefix);
  strcat(fname, localULIgraph);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print out statistics
  maxULI = data.maxULI0;
  maxD = data.maxULI0D;
  if (data.maxULI1 > maxULI) {
    maxULI = data.maxULI1;
    maxD = data.maxULI1D;
  }
  if (data.maxULI5 > maxULI) {
    maxULI = data.maxULI5;
    maxD = data.maxULI5D;
  }
  fprintf(fp, "%lld & %.5f \\\\ \n", to<long long>(maxD), maxULI);
  fclose(fp);

  // ***********************************
  // *** 2. divisors of h statistics ***
  // ***********************************

  // divisors of h for p = 3
  strcpy(fname, fprefix);
  strcat(fname, pdivs3graph);
  if ((fp3 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // divisors of h for p = 5
  strcpy(fname, fprefix);
  strcat(fname, pdivs5graph);
  if ((fp5 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  // divisors of h for p = 7
  strcpy(fname, fprefix);
  strcat(fname, pdivs7graph);
  if ((fp7 = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // calculate pr(p|h)
  sp.reset(3);
  p = sp.next();
  for (i = 0; i < MAXP; i++) {
    prob[i] = nu(p);
    prob[i] = 1 - prob[i];
    p = sp.next();
  }

  // print points to graph
  fprintf(fp3, "%lld\t%.5f\n", to<long long>(gthresh),
          to<double>(to<RR>(data.pdivs0[1][0] + data.pdivs1[1][0]) /
                     (to<RR>(data.total0 + data.total1) * prob[0])));
  fprintf(fp5, "%lld\t%.5f\n", to<long long>(gthresh),
          to<double>(to<RR>(data.pdivs0[2][0] + data.pdivs1[2][0]) /
                     (to<RR>(data.total0 + data.total1) * prob[1])));
  fprintf(fp7, "%lld\t%.5f\n", to<long long>(gthresh),
          to<double>(to<RR>(data.pdivs0[3][0] + data.pdivs1[3][0]) /
                     (to<RR>(data.total0 + data.total1) * prob[2])));

  // close files
  fclose(fp3);
  fclose(fp5);
  fclose(fp7);

  // **********************************************
  // *** 3a. noncyclic odd parts of class group ***
  // **********************************************

  // noncyclic points
  strcpy(fname, fprefix);
  strcat(fname, ncycgraph);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  // print points to graph
  fprintf(fp, "%lld\t%.5f\n", to<long long>(gthresh),
          to<double>((1 - (to<RR>(data.noncyc0 + data.noncyc1) /
                           to<RR>(data.total0 + data.total1))) /
                     cycprob));

  // close file
  fclose(fp);

  // ***********************
  // *** 3b. p-rank data ***
  // ***********************

  // check for ranks 2, 3
  for (i = 2; i <= 4; i++) {

    if (i == 2) {
      // p = 2, r = 2
      strcpy(fname, fprefix);
      strcat(fname, prank22graph);
      if ((fp = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }

      // p = 3, r = 2
      strcpy(fname, fprefix);
      strcat(fname, prank32graph);
      if ((fp3 = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }

      // p = 5, r = 2
      strcpy(fname, fprefix);
      strcat(fname, prank52graph);
      if ((fp5 = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }

      // p = 7, r = 2
      strcpy(fname, fprefix);
      strcat(fname, prank72graph);
      if ((fp7 = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }
    } else if (i == 3) {
      // p = 2, r = 3
      strcpy(fname, fprefix);
      strcat(fname, prank23graph);
      if ((fp = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }

      // p = 3, r = 3
      strcpy(fname, fprefix);
      strcat(fname, prank33graph);
      if ((fp3 = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }

      // p = 5, r = 3
      strcpy(fname, fprefix);
      strcat(fname, prank53graph);
      if ((fp5 = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }

      // p = 7, r = 3
      strcpy(fname, fprefix);
      strcat(fname, prank73graph);
      if ((fp7 = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }
    } else {
      // p = 2, r = 3
      strcpy(fname, fprefix);
      strcat(fname, prank24graph);
      if ((fp = fopen(fname, "a")) == NULL) {
        printf("Cannot open %s.\n", fname);
        exit(1);
      }
    }

    sp.reset(2);
    for (j = 0; j < 4; j++) {
      p = sp.next();

      // calculate ratios
      ratio[j] = nu(p);
      ratio[j] /= pow(to<double>(p), i * i);
      value = nu(p, i);
      ratio[j] /= (value * value);

      // find exponents
      pr0 = pr1 = 0;
      if (i == 2) {
        s3 = s4 = s5 = s6 = 0;
        f3 = f4 = f5 = f6 = 1;
      } else if (i == 3) {
        f3 = MAXE3;
        s3 = f4 = f5 = f6 = 1;
        s4 = s5 = s6 = 0;
      } else if (i == 4) {
        f3 = MAXE3;
        f4 = MAXE4;
        s3 = s4 = f5 = f6 = 1;
        s5 = s6 = 0;
      } else if (i == 5) {
        f3 = MAXE3;
        f4 = MAXE4;
        f5 = MAXE5;
        s3 = s4 = s5 = f6 = 1;
        s6 = 0;
      } else {
        f3 = MAXE3;
        f4 = MAXE4;
        f5 = MAXE5;
        f6 = MAXE6;
        s3 = s4 = s5 = s6 = 1;
      }

      for (e1 = 1; e1 < MAXE1; e1++)
        for (e2 = 1; e2 < MAXE2; e2++)
          for (e3 = s3; e3 < f3; e3++)
            for (e4 = s4; e4 < f4; e4++)
              for (e5 = s5; e5 < f5; e5++)
                for (e6 = s6; e6 < f6; e6++) {
                  pr0 += data.prank0[j][e1][e2][e3][e4][e5][e6];
                  pr1 += data.prank1[j][e1][e2][e3][e4][e5][e6];
                }

      // print points to graph
      if (j == 0)
        fprintf(fp, "%lld\t%.5f\n", to<long long>(gthresh),
                to<double>(to<RR>(pr0 + pr1) /
                           (to<RR>(data.total0 + data.total1) * ratio[0])));
      if (i < 4) {
        if (j == 1)
          fprintf(fp3, "%lld\t%.5f\n", to<long long>(gthresh),
                  to<double>(to<RR>(pr0 + pr1) /
                             (to<RR>(data.total0 + data.total1) * ratio[1])));
        else if (j == 2)
          fprintf(fp5, "%lld\t%.5f\n", to<long long>(gthresh),
                  to<double>(to<RR>(pr0 + pr1) /
                             (to<RR>(data.total0 + data.total1) * ratio[2])));
        else
          fprintf(fp7, "%lld\t%.5f\n", to<long long>(gthresh),
                  to<double>(to<RR>(pr0 + pr1) /
                             (to<RR>(data.total0 + data.total1) * ratio[3])));
      }
    }
    fclose(fp);
    if (i < 4) {
      fclose(fp3);
      fclose(fp5);
      fclose(fp7);
    }
  }

  strcpy(fname, fprefix);
  strcat(fname, maxpgraph);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "%lld\t%lld\t%ld\n", to<long long>(gthresh),
          to<long long>(data.localmaxpD), data.localmaxp);
  fclose(fp);

  return 0;
}

// function to finish and close latex files
void close_output_files() {

  char fname[100];
  FILE *fp;
  cout << "Closing output files ..." << endl;

  // ****************************
  // *** 1. L(1,X) statistics ***
  // ****************************

  // successive L(1) minima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, minL0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $L(1,\\chi)$ minima, $\\Delta=0$ (mod 4). "
              "\\label{table:L0min}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive L(1) minima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minL1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $L(1,\\chi)$ minima, $\\Delta=1$ (mod 8). "
              "\\label{table:L1min}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive L(1) minima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minL5file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $L(1,\\chi)$ minima, $\\Delta=5$ (mod 8). "
              "\\label{table:L5min}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive L(1) maxima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, maxL0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $L(1,\\chi)$ maxima, $\\Delta=0$ (mod 4). "
              "\\label{table:L0max}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive L(1) maxima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxL1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $L(1,\\chi)$ maxima, $\\Delta=1$ (mod 8). "
              "\\label{table:L1max}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive L(1) maxima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxL5file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $L(1,\\chi)$ maxima, $\\Delta=5$ (mod 8). "
              "\\label{table:L5max}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive LLI minima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, minLLI0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $LLI$ minima, $\\Delta=0$ (mod 4). "
              "\\label{table:LLI0min}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive LLI minima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minLLI1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $LLI$ minima, $\\Delta=1$ (mod 8). "
              "\\label{table:LLI1min}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive LLI minima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, minLLI5file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $LLI$ minima, $\\Delta=5$ (mod 8). "
              "\\label{table:LLI5min}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive ULI maxima (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, maxULI0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $ULI$ maxima, $\\Delta=0$ (mod 4). "
              "\\label{table:ULI0max}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive ULI maxima (D = 1 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxULI1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $ULI$ maxima, $\\Delta=1$ (mod 8). "
              "\\label{table:ULI1max}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // successive ULI maxima (D = 5 (mod 8))
  strcpy(fname, fprefix);
  strcat(fname, maxULI5file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }

  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Successive $ULI$ maxima, $\\Delta=5$ (mod 8). "
              "\\label{table:ULI5max}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // local LLI minima
  strcpy(fname, fprefix);
  strcat(fname, localLLIgraph);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fclose(fp);

  // local ULI maxima
  strcpy(fname, fprefix);
  strcat(fname, localULIgraph);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fclose(fp);

  // file of average L function
  strcpy(fname, fprefix);
  strcat(fname, aveLfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Average $L(1,\\chi)$. \\label{table:aveL}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // ***********************************
  // *** 2. divisors of h statistics ***
  // ***********************************

  // file of divisors of h statistics
  strcpy(fname, fprefix);
  strcat(fname, pdivsfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp,
          "\\caption{Number of $h$ divisible by p. \\label{table:pdivs}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h statistics (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivs0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of $h$ divisible by p, $D=0$ (mod 4). "
              "\\label{table:pdivs0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h statistics (D = 1 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivs1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of $h$ divisible by p, $D=1$ (mod 4). "
              "\\label{table:pdivs1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h percentages
  strcpy(fname, fprefix);
  strcat(fname, pdivsperfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(
      fp,
      "\\caption{Percentage of $h$ divisible by p. \\label{table:pdivsper}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h percentages (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivsper0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Percentage of $h$ divisible by p, $D=0$ (mod 4). "
              "\\label{table:pdivsper0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h percentages (D = 1 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivsper1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Percentage of $h$ divisible by p, $D=1$ (mod 4). "
              "\\label{table:pdivsper1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivsratfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of p_p(x). \\label{table:pdivsrat}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h ratios (D = 0 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivsrat0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(
      fp,
      "\\caption{Values of p_p(x), $D=0$ (mod 4). \\label{table:pdivsrat0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of divisors of h ratios (D = 1 (mod 4))
  strcpy(fname, fprefix);
  strcat(fname, pdivsrat1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(
      fp,
      "\\caption{Values of p_p(x), $D=1$ (mod 4). \\label{table:pdivsrat1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of square divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivssqrfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of $p_{p^2}(x)$ \\label{table:pdivssqr0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of square divisors of h ratios, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqr0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of p_{p^2}(x), $D=0$ (mod 4). "
              "\\label{table:pdivssqr0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of square divisors of h ratios, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivssqr1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of p_{p^2}(x), $D=1$ (mod 4). "
              "\\label{table:pdivssqr1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of cube divisors of h ratios
  strcpy(fname, fprefix);
  strcat(fname, pdivscubfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of $p_{p^3}(x)$\\label{table:pdivscub0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of cube divisors of h ratios, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivscub0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of p_{p^3}(x), $D=0$ (mod 4). "
              "\\label{table:pdivscub0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of cube divisors of h ratios, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, pdivscub1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of p_{p^3}(x), $D=1$ (mod 4). "
              "\\label{table:pdivscub1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // *********************************************
  // *** 3a. Noncyclic odd part of class group ***
  // *********************************************

  // file of maximum class numbers
  strcpy(fname, fprefix);
  strcat(fname, maxhfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Maximum class numbers. \\label{table:ncyc}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of maximum odd parts of class numbers
  strcpy(fname, fprefix);
  strcat(fname, maxhoddfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Maximum odd parts of class numbers. "
              "\\label{table:maxhodd}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of noncyclic odd part of class group
  strcpy(fname, fprefix);
  strcat(fname, ncycfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of noncyclic odd parts of class groups. "
              "\\label{table:ncyc}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of noncyclic odd part of class group, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, ncyc0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of noncyclic odd parts of class groups, "
              "$\\Delta=0$ (mod 4). \\label{table:ncyc0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of noncyclic odd part of class group, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, ncyc1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of noncyclic odd parts of class groups, "
              "$\\Delta=1$ (mod 4). \\label{table:ncyc1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // *****************************
  // *** 3b. p-rank statistics ***
  // *****************************

  // file of p-rank statistics
  strcpy(fname, fprefix);
  strcat(fname, prankfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of discriminants with $p$-rank$=r$. "
              "\\label{table:prank}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics, D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prank0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of discriminants with $p$-rank$=r$, "
              "$\\Delta=0$ (mod 4). \\label{table:prank0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics, D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prank1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Number of discriminants with $p$-rank$=r$, "
              "$\\Delta=1$ (mod 4). \\label{table:prank1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics (percentages)
  strcpy(fname, fprefix);
  strcat(fname, prankperfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Percentage of discriminants with $p$-rank$=r$. "
              "\\label{table:prankper}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics (percentages), D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankper0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Percentage of discriminants with $p$-rank$=r$, "
              "$\\Delta=0$ (mod 4). \\label{table:prankper0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics (percentages), D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankper1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Percentage of discriminants with $p$-rank$=r$, "
              "$\\Delta=1$ (mod 4). \\label{table:prankper1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics (ratios)
  strcpy(fname, fprefix);
  strcat(fname, prankratfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of $p_{p,r}(x)$. \\label{table:prankrat}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics (ratios), D = 0 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankrat0file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of $p_{p,r}(x)$, $\\Delta=0$ (mod 4). "
              "\\label{table:prankrat0}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of p-rank statistics (ratios), D = 1 (mod 4)
  strcpy(fname, fprefix);
  strcat(fname, prankrat1file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of $p_{p,r}(x)$, $\\Delta=1$ (mod 4). "
              "\\label{table:prankrat1}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // ***************************
  // *** 4. First occurences ***
  // ***************************

  // files of first occurences of rank 2 2-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, first22file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 2 2-Sylow subgroups. "
              "\\label{table:first22}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 2 p-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, firstp2file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 2 $p$-Sylow subgroups. "
              "\\label{table:firstp2}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 3 2-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, first23file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 3 2-Sylow subgroups. "
              "\\label{table:first23}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 3 p-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, firstp3file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 3 $p$-Sylow subgroups. "
              "\\label{table:firstp3}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 4 2-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, first24file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 4 2-Sylow subgroups. "
              "\\label{table:first24}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 4 p-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, firstp4file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 4 $p$-Sylow subgroups. "
              "\\label{table:firstp4}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 5 2-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, first25file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 5 2-Sylow subgroups. "
              "\\label{table:first25}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 5 p-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, firstp5file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 5 $p$-Sylow subgroups. "
              "\\label{table:firstp5}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 6 2-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, first26file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 6 2-Sylow subgroups. "
              "\\label{table:first26}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences of rank 6 p-Sylow groups
  strcpy(fname, fprefix);
  strcat(fname, firstp6file);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Non-cyclic rank 6 $p$-Sylow subgroups. "
              "\\label{table:firstp6}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences doubly noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, twononcycfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(
      fp,
      "\\caption{Doubly non-cyclic class groups. \\label{table:twononcyc}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences trebly noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, threenoncycfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Trebly non-cyclic class groups. "
              "\\label{table:threenoncyc}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // files of first occurences quadruply noncyclic groups
  strcpy(fname, fprefix);
  strcat(fname, fournoncycfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Quadruply non-cyclic class groups. "
              "\\label{table:fournoncyc}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // ****************************************
  // *** 5. number of generators required ***
  // ****************************************

  // file of maximum prime ideal required
  strcpy(fname, fprefix);
  strcat(fname, maxpfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Maximum prime ideal required. \\label{table:maxp}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of maximum prime ideal ratios
  strcpy(fname, fprefix);
  strcat(fname, maxratfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{Values of $\\max{p}/\\log{\\Delta}$ and "
              "$\\max{p}/\\log^2{\\Delta}$. \\label{table:maxrat}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of psplit
  strcpy(fname, fprefix);
  strcat(fname, psplitfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{First $\\Delta$ needing odd prime splits up to $p$. "
              "\\label{table:psplit}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);

  // file of kpsplit
  strcpy(fname, fprefix);
  strcat(fname, kpsplitfile);
  if ((fp = fopen(fname, "a")) == NULL) {
    printf("Cannot open %s.\n", fname);
    exit(1);
  }
  fprintf(fp, "\\end{tabular}\n");
  fprintf(fp, "\\end{center}\n");
  fprintf(fp, "\\caption{First $\\Delta$ needing $k$ odd prime splits. "
              "\\label{table:kpsplit}}\n");
  fprintf(fp, "\\end{table}\n");
  fclose(fp);
}

double nu(long p) {

  RR tolerance = power(to<RR>(0.1), 20);
  RR diff, product, oldproduct, term;
  // double product, oldproduct, term;
  double pk, pr;

  set(product);
  set(diff);
  set(pk);
  pr = to<double>(p);

  while (diff > tolerance) {
    pk *= pr;
    term = to<RR>(1.0 / pk);
    term = 1.0 - term;

    product *= term;
    diff = abs(product - oldproduct);
    oldproduct = product;
  }
  return to<double>(product);
}

double nu(long p, long alpha) {

  long k;
  RR product, term;
  RR pk, pr;

  set(product);
  set(pk);
  pr = to<RR>(p);
  for (k = 1; k <= alpha; k++) {
    pk *= pr;
    term = 1.0 / pk;
    term = 1.0 - term;
    product *= term;
  }
  return to<double>(product);
}
