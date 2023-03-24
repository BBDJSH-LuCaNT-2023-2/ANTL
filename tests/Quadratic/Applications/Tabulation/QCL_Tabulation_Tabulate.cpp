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

#include "timer.hpp"

#include <fstream>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>

#ifndef INVARIANTS_TABULATION_TEST
#define INVARIANTS_TABULATION_TEST

// #define LIST_SIZE_QUADRATIC 33554432
#define LIST_SIZE_QUADRATIC 1000000

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

//     long x;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (myrank == 0) {
      ZZ total_n, total_t1, total_tr, total_tc, IS;
      long type, L, H, i, idx, n, t1, tr, tc, jobs, saveH, savenf;
      bool done;
      char fname[50], hname[50], prefix[100]; //,mkdir[50],check[50];
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

      size_t len = sizeof(char)*50;
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

      for (i = L; i < H; ++i) {
        finished[i] = false;
      }

      sprintf(fname, "check-%ld-%s", L, argv[3]);
      cout << "checkpoint file = " << fname << endl << endl;

      infile.open(fname);
      if (infile.fail()) {
        clear(total_n);
        clear(total_t1);
        clear(total_tr);
        clear(total_tc);
      } else {
        infile >> savenf;
        infile >> saveH;
        if (H >= saveH) {
          infile >> total_n;
          infile >> total_t1;
          infile >> total_tr;
          infile >> total_tc;
          for (i = 0; i < savenf; ++i) {
            infile >> finished[L + i];
          }
        } else {
          clear(total_n);
          clear(total_t1);
          clear(total_tr);
          clear(total_tc);
        }
        infile.close();
      }

      saveH = H;
      savenf = numfinished;

      // continue to collect stats and start new jobs until all are finished
      done = true;
      for (i = L; i < H && done; ++i) {
        done = finished[i];
      }
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
          while (idx < H && finished[idx]) {
            ++idx;
          }
        }

        while (idx < H) {
          // receive computational data, start new job if necessary
          MPI_Recv(buffer, sizeof(char)*BUFLEN, MPI_PACKED, MPI_ANY_SOURCE, TIME_DATA,
                   MPI_COMM_WORLD, &rc);
          --jobs;

          position = 0;
          MPI_Unpack(buffer, BUFLEN, &position, &i, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &n, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &t1, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &tr, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &tc, 1, MPI_LONG,
                     MPI_COMM_WORLD);

          cout << "Interval " << i << " - " << n << " fields, Total Time: " << flush;
          MyTime(t1*10000);
          cout << " (Regulator: " << flush; MyTime(tr);
          cout << " Class Group: " << flush; MyTime(tc);
          cout << ")" << endl;

          total_n += n;
          total_t1 += t1;
          total_tr += tr;
          total_tc += tc;

          // update checkpoint file
          finished[i] = true;
          outfile.open(fname);
          outfile << savenf << " " << saveH << " ";
          outfile << total_n << " " << total_t1 << " " << total_tr << " " << total_tc << " ";
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
          MPI_Unpack(buffer, BUFLEN, &position, &tr, 1, MPI_LONG,
                     MPI_COMM_WORLD);
          MPI_Unpack(buffer, BUFLEN, &position, &tc, 1, MPI_LONG,
                     MPI_COMM_WORLD);

          cout << "Interval " << i << " - " << n << " fields, Total Time: " << flush;
          MyTime(t1*10000);
          cout << " (Regulator Time: " << flush; MyTime(tr);
          cout << " Class Group Time: " << flush; MyTime(tc);
          cout << endl;

          total_n += n;
          total_t1 += t1;
          total_tr += tr;
          total_tc += tc;

          // update checkpoint file
          finished[i] = true;
          outfile.open(fname);
          outfile << savenf << " " << saveH << " ";
          outfile << total_n << " " << total_t1 << " " << total_tr << " " << total_tc << " ";
          for (i = L; i < H; ++i)
            outfile << finished[i] << " ";
          outfile << endl;
          outfile.close();
        }

        done = true;
        for (i = L; i < H && done; ++i)
          done = finished[i];
      }

      std::cout << "send a \"finished\" flag to all slaves" << std::endl;
      // send a "finished" flag to all slaves
      idx = -1;
      for (i = 1; i <= machines; ++i)
        MPI_Send(&idx, 1, MPI_LONG, i, INTERVAL_DATA, MPI_COMM_WORLD);

      // output runtime statistics
      cout << "\n\nSummary:" << endl;
      cout << "Total Fields:  " << total_n << endl;
      cout << "Total time:  " << flush;
      MyTime(total_t1*10000);
      cout << endl;
      cout << "Total regulator time:  " << flush;
      MyTime(total_tr);
      cout << endl;
      cout << "Total classgroup time:  " << flush;
      MyTime(total_tc);
      cout << endl;
      cout << "Real time:  " << flush;
      MyTime((total_t1*10000) / machines);
      cout << endl;
      cout << "Real regulator time:  " << flush;
      MyTime(total_tr / machines);
      cout << endl;
      cout << "Real class group time:  " << flush;
      MyTime(total_tc / machines);
      cout << endl;
      cout << "Avg per field:  " << flush;
      MyTime(total_t1 / total_n);
      cout << endl << endl;
    } else {
      ZZ Dlist[LIST_SIZE_QUADRATIC], L, H, maxH, IS;
//       QuadraticOrder<ZZ> QO;
//       vec_ZZ Cl;
      long n, i, j, idx, max_idx, Dl, oldDl, t1, tc = 0, rank;
      vector<long> tr;
      long long D, oldD;
      char fname[50], hname[50], zipper[100], check[100], mkdir[100];
      long alg = 1;
      timer t;
      ofstream outfile;
      MPI_Status rc;
      char buffer[BUFLEN];
      int position;

      //
      // Initializations
      //
      size_t len = sizeof(char)*50;
      gethostname(hname, len);

      cout << "Proc " << myrank << " (" << hname << "):  Initializing" << endl;

      // turn off verbosity for h computations
      // quadratic_order<long>::verbose(0);

      long type = atol(argv[1]);
      conv(IS, argv[4]);
      maxH = to<ZZ>(atol(argv[3])) * IS;

      cout << "Proc " << myrank << " (" << hname << "):  Initializing prime list" << endl;
      // initialize prime list for discriminant generation
      if (type == 0) {
//         get_Dlist_imag(maxH, maxH, Dlist, n, 1);
      }
      else {
        get_Dlist_real(maxH, maxH, Dlist, n, 1);
      }

      // initialize L(1,X) tables and prime list for discriminant generation
//       QO.use_Lfunction_tables(to<long>(maxH));

//       cout << "Proc " << myrank << " (" << hname << "): check directory" << endl;
      // check directory
      struct stat buf;
      sprintf(check, "%s/%s", argv[5], argv[6]);
      if (stat(check, &buf)) {
        sprintf(mkdir, "mkdir %s/%s", argv[5], argv[6]);
        system(mkdir);
      }

      cout << "Proc " << myrank << " (" << hname << "): main loop - process interval idx until idx < 0 received" << endl;
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

        if (type == 0) {
          sprintf(fname, "%s/%s/iq-data.%ld", argv[5], argv[6], idx);
        }
        else {
          sprintf(fname, "%s/%s/rq-data.%ld", argv[5], argv[6], idx);
        }
        outfile.open(fname);

        t.start_timer();
        tr = {0, 0, 0, 0, 0};
        tc = 0;

        // compute Cl for all fundamental discriminants D with L <= |D| <= H
        // File format:
        //   L H n (n = number of fields)
        //   D #primes pmax h rank C[0] Cl[1] ... CL[rank-1]
        //   ...
        //   time (in csec)

        if (type == 0) {
//           get_Dlist_imag(L, H, Dlist, n);
        }
        else {
          if (L == 1) {
//             get_Dlist_real(to_ZZ(2), H, Dlist, n);
            get_Dlist_real(to_ZZ(2), H, Dlist, n, 0);
          }
          else {
//             get_Dlist_real(L, H, Dlist, n);
            get_Dlist_real(L, H, Dlist, n, 0);
          }
        }

        outfile << L << " " << H << " " << n << endl;

        // set global # of terms and prec for L-function approximations for this
        // interval
//         QO.set_Lfunction_global(to<long>(H));

        oldDl = to<long>(L);
        for (i = 0; i < n; ++i) {
//           std::cout << "Dlist[i] is " << Dlist[i] << std::endl;
          conv(Dl, Dlist[i]);
//           QO.assign(Dl);
          QuadraticOrder<long> QO{Dl};
          outfile << Dl << flush;

          // ===TEMPORARY INITIALIZATION===
          // TODO all of this initialization should be done from within QuadraticOrder
          // using new
          QuadraticNumber<long> quad_number1{QO};
          QuadraticNumber<long> quad_number2{QO};
          QuadraticNumber<long> quad_number3{QO};

          MultiplyNucompOpt<long> mul_nucomp_opt_object{QO};
          mul_nucomp_opt_object.set_RelativeGenerator(quad_number1);
          QO.set_mul_nucomp_opt(mul_nucomp_opt_object);

          SquareNuduplOpt<long> sqr_nudupl_opt_object{QO};
          sqr_nudupl_opt_object.set_RelativeGenerator(quad_number2);
          QO.set_sqr_best(sqr_nudupl_opt_object);

          ReducePlainRealOpt<long> red_plain_real_opt_object{};
          red_plain_real_opt_object.set_RelativeGenerator(quad_number3);
          QO.set_red_best(red_plain_real_opt_object);

          L_function<long> l_function;
          l_function.init(to_long(Dl), 2);

          tuple<double, ZZ, vector<long>> regulator_and_hstar_tuple = get_regulator_and_hstar(QO, l_function);

          double regulator = std::get<0>(regulator_and_hstar_tuple);
          ZZ h_star = std::get<1>(regulator_and_hstar_tuple);
          for(int i = 0; i < 5; i++) {
            tr[i] += std::get<2>(regulator_and_hstar_tuple)[i];
          }

          vector<long> class_group;
          if(alg == 0) {
            // Computing and timing a single class group BSGS computation
            std::pair<vector<long>, long> class_group_pair = get_class_group_BSGS(QO, regulator, h_star);
            class_group = std::get<0>(class_group_pair);
            tc += std::get<1>(class_group_pair);
          }
          else if(alg == 1) {
            // Computing and timing a single class group BS computation
            std::pair<vector<long>, long> class_group_pair = get_class_group_BS(QO, regulator, h_star);
            class_group = std::get<0>(class_group_pair);
            tc += std::get<1>(class_group_pair);
          }

          //Formatting the class group first
          class_group.erase(std::remove(class_group.begin(), class_group.end(), 1), class_group.end());

          //Outputting to stream
//           std::cout << " " << regulator*1000 << flush;
          outfile << std::fixed;
          outfile << " " << regulator << flush;
          outfile << " " << class_group << std::endl;
        }

        t.stop_timer();
        t1 = t.user_time();

//         outfile << t1 << endl;
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
        MPI_Pack(&tr[0], 1, MPI_LONG, buffer, BUFLEN, &position, MPI_COMM_WORLD);
        MPI_Pack(&tc, 1, MPI_LONG, buffer, BUFLEN, &position, MPI_COMM_WORLD);
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
