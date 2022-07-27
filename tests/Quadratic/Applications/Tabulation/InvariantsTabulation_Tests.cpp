/**
 * @file qcl_table.cpp
 * @author Michael Jacobson
 * @version $Header:
 * /scrinium/ANTL/ANTL/appl/quadratic/iqcl-table.cpp,v
 * 1.2 2013/02/26 06:29:09 jacobs Exp $
 * @brief Quadratic field class group tabulator
 * program.
 *
 * MPI-based client/server program for tabulating
 * all class groups for quadratic fields with
 * discriminant D such that L <= |D| <= H.
 *
 * Usage:  iqcl-table L H IS"
 *   L - lower index for interval
 *   H - upper index for interval
 *   IS - interval width
 *   Computes Cl for all D with L*IS <= D <= H*IS
 *   Must have IS <= INT_SIZE_QUADRATIC
 */

#define INTERVAL_DATA 0
#define TIME_DATA 1
#define BUFLEN 50

// #include <ANTL/quadratic/quadratic_order.hpp>
#include <fstream>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef INVARIANTS_TABULATION_TEST
#define INVARIANTS_TABULATION_TEST

#define LIST_SIZE_QUADRATIC 33554432

#include "AuxillaryFunctions.hpp"

#include <mpi.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <string>

#include <ANTL/Quadratic/ClassGroup/ClassGroupBSGSReal.hpp>
#include <ANTL/Quadratic/Regulator/RegulatorLenstra_ZZ.hpp>

using namespace NTL;
using namespace ANTL;
// 1 thread: 58961 ms
// 12 thread: 243066 ms

  int main(int argc, char **argv) {
    int myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) {
      ZZ total_n, total_t1, IS;
      register long type, L, H, i, idx, n, t1, jobs, saveH, savenf;
      bool done;
      char fname[50], hname[50], prefix[50]; //,mkdir[50],check[50];
      ifstream infile;
      ofstream outfile;

      char buffer[BUFLEN];
      int temp, position = 0;
      MPI_Status rc;

      // get number of available machines
      MPI_Comm_size(MPI_COMM_WORLD, &temp);
      long machines = (long)temp;

      --machines;

      if (argc != 7) {
        cerr << "Usage:  qcl-table type L H IS prefix dir" << endl;
        cerr << " type - 0 (imaginary) or 1 (real)" << endl;
        cerr << " L - lower index for interval" << endl;
        cerr << " H - upper index for interval" << endl;
        cerr << " IS - interval width" << endl;
        cerr << " prefix - prefix for nodes' local data files" << endl;
        cerr << " dir - directory for data files" << endl;
        cerr << " Computes Cl for all fundamental D with L*IS <= D <= H*IS"
             << endl;
        cerr << " Must have IS <= LIST_SIZE_QUADRATIC = " << LIST_SIZE_QUADRATIC
             << endl;
        exit(1);
      }

      // get range of interval indices
      type = atol(argv[1]);
      L = atol(argv[2]);
      H = atol(argv[3]);
      conv(IS, argv[4]);

      size_t len = 50;
      gethostname(hname, len);

      cout << "Master:  " << hname << endl;
      cout << "Computing CL for all fundamental D in " << to<ZZ>(L) * IS
           << " to " << to<ZZ>(H) * IS << endl;
      if (type == 0)
        cout << "tablating imaginary fields" << endl;
      else
        cout << "tablating real fields" << endl;
      cout << "indices " << L << " to " << H << endl;
      cout << "interval width:  " << IS << endl;
      cout << "local file prefix:  " << argv[5] << endl;
      cout << "directory: " << argv[6] << endl;
      cout << "Using " << machines << " CPUs" << endl;
      cout << endl;
      cout << flush;

      // read checkpoint file, if it exists
      bool finished[1 + H];
      long numfinished = H - L;

      for (i = L; i < H; ++i)
        finished[i] = false;

      sprintf(fname, "check-%ld-%s", L, argv[3]);
      cout << "checkpoint file = " << fname << endl << endl;

      infile.open(fname);
      if (infile.fail()) {
        clear(total_n);
        clear(total_t1);
      } else {
        infile >> savenf;
        infile >> saveH;
        if (H >= saveH) {
          infile >> total_n;
          infile >> total_t1;
          for (i = 0; i < savenf; ++i)
            infile >> finished[L + i];
        } else {
          clear(total_n);
          clear(total_t1);
        }
        infile.close();
      }

      saveH = H;
      savenf = numfinished;

      // continue to collect stats and start new jobs until all are finished
      done = true;
      for (i = L; i < H && done; ++i)
        done = finished[i];

      while (!done) {
        // start first jobs on all machines
        idx = L;
        while (idx < H && finished[idx])
          ++idx;

        jobs = 0;
        for (i = 1; i <= machines && idx < H; ++i) {
          // get next interval to compute (idx < 0 indicates termination)
          MPI_Send(&idx, 1, MPI_LONG, i, INTERVAL_DATA, MPI_COMM_WORLD);
          ++jobs;

          ++idx;
          while (idx < H && finished[idx])
            ++idx;
        }

        while (idx < H) {
          // receive computational data, start new job if necessary
          MPI_Recv(buffer, BUFLEN, MPI_PACKED, MPI_ANY_SOURCE, TIME_DATA,
                   MPI_COMM_WORLD, &rc);
          --jobs;

          position = 0;
          MPI_Unpack(buffer, BUFLEN, &position, &i, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &n, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &t1, 1, MPI_LONG,
                     MPI_COMM_WORLD);

          cout << "Interval " << i << " - " << n << " fields, time " << flush;
          MyTime(t1);
          cout << endl;

          total_n += n;
          total_t1 += t1;

          // update checkpoint file
          finished[i] = true;
          outfile.open(fname);
          outfile << savenf << " " << saveH << " ";
          outfile << total_n << " " << total_t1 << " ";
          for (i = L; i < H; ++i)
            outfile << finished[i] << " ";
          outfile << endl;
          outfile.close();

          // send next interval to compute
          while (idx < H && finished[idx])
            ++idx;
          if (idx < H) {
            MPI_Send(&idx, 1, MPI_LONG, rc.MPI_SOURCE, INTERVAL_DATA,
                     MPI_COMM_WORLD);
            ++jobs;
            ++idx;
          }
        }

        // get final results
        while (jobs > 0) {
          // receive computational data, start new job if necessary
          MPI_Recv(buffer, BUFLEN, MPI_PACKED, MPI_ANY_SOURCE, TIME_DATA,
                   MPI_COMM_WORLD, &rc);
          --jobs;

          position = 0;
          MPI_Unpack(buffer, BUFLEN, &position, &i, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &n, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &t1, 1, MPI_LONG,
                     MPI_COMM_WORLD);

          cout << "Interval " << i << " - " << n << " fields, time " << flush;
          MyTime(t1);
          cout << endl;

          total_n += n;
          total_t1 += t1;

          // update checkpoint file
          finished[i] = true;
          outfile.open(fname);
          outfile << savenf << " " << saveH << " ";
          outfile << total_n << " " << total_t1 << " ";
          for (i = L; i < H; ++i)
            outfile << finished[i] << " ";
          outfile << endl;
          outfile.close();
        }

        done = true;
        for (i = L; i < H && done; ++i)
          done = finished[i];
      }

      // send a "finished" flag to all slaves
      idx = -1;
      for (i = 1; i <= machines; ++i)
        MPI_Send(&idx, 1, MPI_LONG, i, INTERVAL_DATA, MPI_COMM_WORLD);

      // output runtime statistics
      cout << "\n\nSummary:" << endl;
      cout << "Total Fields:  " << total_n << endl;
      cout << "Total time:  " << flush;
      MyTime(total_t1);
      cout << endl;
      cout << "Real time:  " << flush;
      MyTime(total_t1 / machines);
      cout << endl;
      cout << "Avg per field:  " << flush;
      MyTime(total_t1 / total_n);
      cout << endl << endl;
    } else {

      ZZ Dlist[LIST_SIZE_QUADRATIC], L, H, maxH, IS;
      quadratic_order<long> QO;
      vec_ZZ Cl;
      register long n, i, j, idx, max_idx, Dl, oldDl, t1, rank;
      long long D, oldD;
      char fname[50], hname[50], zipper[100], check[50], mkdir[50];
      timer t;
      ofstream outfile;
      MPI_Status rc;
      char buffer[BUFLEN];
      int position;

      //
      // Initializations
      //

      size_t len = 100;
      gethostname(hname, len);

      cout << "Proc " << myrank << " (" << hname << "):  Initializing" << endl;

      // turn off verbosity for h computations
      quadratic_order<long>::verbose(0);

      long type = atol(argv[1]);
      conv(IS, argv[4]);
      maxH = to<ZZ>(atol(argv[3])) * IS;

      // initialize prime list for discriminant generation
      if (type == 0)
        get_Dlist_imag(maxH, maxH, Dlist, n, 1);
      else
        get_Dlist_real(maxH, maxH, Dlist, n, 1);

      // initialize L(1,X) tables and prime list for discriminant generation
      QO.use_Lfunction_tables(to<long>(maxH));

      // check directory
      struct stat buf;
      sprintf(check, "%s/%s", argv[5], argv[6]);
      if (stat(check, &buf)) {
        sprintf(mkdir, "mkdir %s/%s", argv[5], argv[6]);
        system(mkdir);
      }

      //
      // main loop - process interval idx until idx < 0 received
      //

      // get interval to compute (idx < 0 indicates termination)
      MPI_Recv(&idx, 1, MPI_LONG, 0, INTERVAL_DATA, MPI_COMM_WORLD, &rc);

      while (idx >= 0) {
        // get L and H
        L = to_ZZ(idx) * IS;
        H = L + IS;
        ++L;

        cout << "Proc " << myrank << " (" << hname << "):  " << idx << ", " << L
             << " to " << H << endl;

        if (type == 0)
          sprintf(fname, "%s/%s/iq-data.%ld", argv[5], argv[6], idx);
        else
          sprintf(fname, "%s/%s/rq-data.%ld", argv[5], argv[6], idx);
        outfile.open(fname);

        t.start_timer();

        // compute Cl for all fundamental discriminants D with L <= |D| <= H
        // File format:
        //   L H n (n = number of fields)
        //   D #primes pmax h rank C[0] Cl[1] ... CL[rank-1]
        //   ...
        //   time (in csec)

        if (type == 0)
          get_Dlist_imag(L, H, Dlist, n);
        else {
          if (L == 1)
            get_Dlist_real(to_ZZ(2), H, Dlist, n);
          else
            get_Dlist_real(L, H, Dlist, n);
        }

        outfile << L << " " << H << " " << n << endl;

        // set global # of terms and prec for L-function approximations for this
        // interval
        QO.set_Lfunction_global(to<long>(H));

        oldDl = to<long>(L);
        for (i = 0; i < n; ++i) {
          conv(Dl, Dlist[i]);
          QO.assign(Dl);
          Cl = QO.class_group(CLASS_GROUP_BSGS);
          rank = QO.get_rank();
          outfile << (labs(Dl) - oldDl) << " " << flush;
          // outfile << Dl << " " << flush;
          outfile << QO.get_nump() << " " << flush;
          outfile << QO.get_pmax() << " " << flush;
          if (type == 1)
            outfile << QO.regulator().get_log() << " " << flush;
          outfile << rank << flush;
          for (j = 0; j < rank; ++j)
            outfile << " " << Cl[j] << flush;
          outfile << endl;
          oldDl = labs(Dl);
        }

        t.stop_timer();
        t1 = t.user_time();

        outfile << t1 << endl;
        outfile.close();

        // gzip the output file
        if (type == 0)
          sprintf(zipper, "gzip %s/%s/iq-data.%ld", argv[5], argv[6], idx);
        else
          sprintf(zipper, "gzip %s/%s/rq-data.%ld", argv[5], argv[6], idx);
        system(zipper);

        // send computation stats to master
        position = 0;
        MPI_Pack(&idx, 1, MPI_LONG, buffer, BUFLEN, &position, MPI_COMM_WORLD);
        MPI_Pack(&n, 1, MPI_LONG, buffer, BUFLEN, &position, MPI_COMM_WORLD);
        MPI_Pack(&t1, 1, MPI_LONG, buffer, BUFLEN, &position, MPI_COMM_WORLD);
        MPI_Send(buffer, BUFLEN, MPI_PACKED, 0, TIME_DATA, MPI_COMM_WORLD);

        // get next interval to compute (idx < 0 indicates termination)
        MPI_Recv(&idx, 1, MPI_LONG, 0, INTERVAL_DATA, MPI_COMM_WORLD, &rc);
      }
    }

    MPI_Finalize();

    if (myrank == 0)
      cout << "Finished!" << endl;

    exit(0);
  }
#endif
