/**
 * @file Test-interval.cpp
 * @author Michael Jacobson
 * @version $Header: /scrinium/ANTL/ANTL/appl/quadratic/qcl_interval.cpp,v 1.8
 * 2013/02/26 06:29:09 jacobs Exp $
 * @brief Quadratic field interval test program.
 *
 * Computes class groups for all quadratic fields of discriminant D
 *   with L <= |D| <= H
 * Usage:  qcl-interval L H type table cond vb
 *   L - lower bound for |D|
 *   H - upper bound for |D|
 *   type - 0 (imaginary), 1 (real)
 *   table - 0 (do not use L(1,X) tables), 1 (use them)
 *   uncond - 0 (conditional on ERH), 1 (unconditional)
 *   vb - verbosity (for QuadraticOrder class)
 *   kind - bsgs or bs
 *   prec - 0 (automatic), 1 (long), 2 (long long), >2 (ZZ)
 */

// #define LIST_SIZE_QUADRATIC 33554432
#define LIST_SIZE_QUADRATIC 1000000

#include "AuxillaryFunctions.hpp"
#include "timer.hpp"


#include <ANTL/Quadratic/QuadraticOrder.hpp>

NTL_CLIENT
using namespace ANTL;

int main(int argc, char **argv) {
  ZZ Dlist[LIST_SIZE_QUADRATIC], L, H, D, pmax;
  vec_ZZ Cl;
//   QuadraticOrder<ZZ> QO;
//   QuadraticOrder<long> QOl;
//   QuadraticOrder<long long> QOll;
  long long Dll;
  long n, i, Dl;
  timer t;
  int vb, tp, table, uncond, method, prec;

  if (argc != 9) {
    cerr << "Usage:  qcl_interval L H type table uncond vb alg prec" << std::endl;
    cerr << "  L - lower bound for |D|" << std::endl;
    cerr << "  H - upper bound for |D|" << std::endl;
    cerr << "  type = 0 (imag), 1 (real)" << std::endl;
    cerr << "  table - 0 (do not use L(s,X) tables), 1 (use them)" << std::endl;
    cerr << "  uncond - 0 (conditional on ERH), 1 (unconditional)" << std::endl;
    cerr << "  vb - verbosity" << std::endl;
    cerr << "  alg - 0 (bsgs), 1 (bs)" << std::endl;
    cerr << "  prec - 0 (automatic), 1 (long), 2 (long long), >2 (ZZ)" << std::endl;
    exit(1);
  }

  conv(L, argv[1]);
  conv(H, argv[2]);
  tp = atoi(argv[3]);

  table = atoi(argv[4]);
  uncond = atoi(argv[5]);
  vb = atoi(argv[6]);
  method = atoi(argv[7]);
  prec = atoi(argv[8]);

  std::cout << "Computing class groups for interval " << L << " <= |D| <= " << H
       << std::endl;
  if (tp == 0)
    std::cout << "Imaginary fields" << std::endl;
  else
    std::cout << "Real fields" << std::endl;
  if (table)
    std::cout << "Using L(s,X) tables" << std::endl;
  if (uncond)
    std::cout << "Using unconditional L(1,X) approximations" << std::endl;
  if (method == 0)
    std::cout << "Using BSGS" << std::endl;
  else
    std::cout << "Using BS" << std::endl;

  // set quadratic order verbosity level
//   QuadraticOrder<ZZ>::verbose(vb);
//   QuadraticOrder<long>::verbose(vb);
//   QuadraticOrder<long long>::verbose(vb);

  // initialize unconditional Lfunction approximations if uncond is set
//   if (uncond) {
//     QOl.set_unconditional();
//     QOll.set_unconditional();
//     QO.set_unconditional();
//   }

  // time computation from this point on
  t.start_timer();

  // generate list of fundamental D with L <= |D| <= H
  if (tp == 0) {
//     get_Dlist_imag(L, H, Dlist, n, 1); // initialize prime list
//     get_Dlist_imag(L, H, Dlist, n);
  } else {
    get_Dlist_real(L, H, Dlist, n, 1); // initialize prime list
//     get_Dlist_real(L, H, Dlist, n);
    get_Dlist_real(L, H, Dlist, n, 0);
  }
  std::cout << "\n# fields = " << n << std::endl;

/*
  if (prec == 1 || (prec == 0 && H < NTL_SP_BOUND)) {
    // use longs
    std::cout << "Using long" << std::endl;

    // initialize table-driven L-function approximations if table is set
    conv(Dl, H);
//     if (table)
//       QOl.use_Lfunction_tables(Dl);
//     else
//       QOl.set_Lfunction_global(Dl);

    // compute class groups
    for (i = 0; i < n; ++i) {
      conv(Dl, Dlist[i]);
      QOl.assign(Dl);
      std::cout << "D = " << Dl << flush;
      if (method == 0)
        Cl = QOl.class_group(CLASS_GROUP_BSGS);
      else
        Cl = QOl.class_group(CLASS_GROUP_BS);

      std::cout << ", # primes = " << QOl.get_nump() << flush; // prime ideals
      std::cout << ", pmax = " << QOl.get_pmax() << flush;
      if (tp > 0)
        std::cout << ", R = " << QOl.regulator().get_log() << flush;
      std::cout << ", h = " << QOl.class_number() << flush;
      std::cout << ", Cl = " << Cl << std::endl;
    }
  } else if (prec == 2 || (prec == 0 && to<long long>(H) <=
                                            to<long long>(NTL_SP_BOUND) *
                                                to<long long>(NTL_SP_BOUND))) {
    // use long longs
    std::cout << "Using long long" << std::endl;

    // initialize table-driven L-function approximations if table is set
    Dll = to<long long>(H);
    if (table)
      QOll.use_Lfunction_tables(Dll);
    else
      QOll.set_Lfunction_global(Dll);

    // compute class groups
    for (i = 0; i < n; ++i) {
      Dll = to<long long>(Dlist[i]);
      QOll.assign(Dll);
      std::cout << "D = " << Dll << flush;
      if (method == 0)
        Cl = QOll.class_group(CLASS_GROUP_BSGS);
      else
        Cl = QOll.class_group(CLASS_GROUP_BS);

      std::cout << ", # primes = " << QOll.get_nump() << flush;
      std::cout << ", pmax = " << QOll.get_pmax() << flush;
      if (tp > 0)
        std::cout << ", R = " << QOll.regulator().get_log() << flush;
      std::cout << ", h = " << QOll.class_number() << flush;
      std::cout << ", Cl = " << Cl << std::endl;
    }
  } else {
*/
    // use ZZ
    std::cout << "Using ZZ" << std::endl;

    // initialize table-driven L-function approximations if table is set
//     if (table)
//       QO.use_Lfunction_tables(H);
//     else
//       QO.set_Lfunction_global(H);

    // compute class groups
    /*
    FILE *out = fopen("sample2.dat","wb");
    long long llL, llH;
    llL = to<long long>(L);
    llH = to<long long>(H);

    fwrite(&llL,sizeof(long long), 1, out);
    fwrite(&llH,sizeof(long long), 1, out);
    fwrite(&n,sizeof(long),1,out);
    */

    std::cout << L << " " << H << " " << n << std::endl;

    ZZ oldD;
    oldD = L;

    for (i = 0; i < n; ++i) {
      D = Dlist[i];
//       QO.assign(D);
      QuadraticOrder<ZZ> QO{ZZ(D)};
      std::cout << "Current discriminant is " << D << std::endl;

      // std::cout << "D = " << D << flush;

//       if (method == 0)
//         Cl = QO.class_group(CLASS_GROUP_BSGS);
//       else
//         Cl = QO.class_group(CLASS_GROUP_BS);
      pair<double, vector<long>> regulator_and_class_group = get_regulator_and_class_group(QO);
      double regulator = regulator_and_class_group.first;
      vector<long> class_group = regulator_and_class_group.second;
      /*
                        unsigned char diff = (unsigned char) to<long>(abs(D) -
         oldD); if ((abs(D) - oldD) >= 256)  std::cout << "ERROR:  D = " << D << ",
         oldD = " << oldD << ", diff = " << abs(D) - oldD << std::endl;
                        fwrite(&diff,sizeof(unsigned char),1,out);

                        unsigned char np = (unsigned char) QO.get_nump();
                        fwrite(&np,sizeof(unsigned char),1,out);

                        short mp = (short) to<long>(QO.get_pmax());
                        fwrite(&mp,sizeof(short),1,out);

                        unsigned char rank = (unsigned char) QO.get_rank();
                        fwrite(&rank,sizeof(unsigned char),1,out);

                        for (long j=0; j<(long) rank; ++j) {
                          long ediv = to<long>(Cl[j]);
                          fwrite(&ediv,sizeof(long),1,out);
                        }

                        quad_float dR =
         to<quad_float>(QO.regulator().get_log());
                        fwrite(&dR,sizeof(quad_float),1,out);

                        oldD = abs(D);
      */

//       long rank = QO.get_rank();
//       std::cout << (abs(D) - oldD) << " " << flush;
//       std::cout << QO.get_nump() << " " << flush;
//       std::cout << QO.get_pmax() << " " << flush;
//       std::cout << rank << flush;
//       for (long j = 0; j < rank; ++j)
//         std::cout << " " << Cl[j] << flush;
      std::cout << " " << class_group << flush;
      std::cout << " " << regulator << flush;
      std::cout << std::endl;
      oldD = abs(D);

      /*
                std::cout << ", # primes = " << QO.get_nump() << flush;
                std::cout << ", pmax = " << QO.get_pmax() << flush;
                if (tp > 0)
                  std::cout << ", R = " << to<quad_float>(QO.regulator().get_log())
         << flush; std::cout << ", h = " << QO.class_number() << flush; std::cout << ", Cl
         = " << Cl << std::endl;
      */
    }
    // fclose(out);
//  }

  // stop timing and output total elapsed CPU time
  t.stop_timer();

  if (vb) {
    std::cout << "\nTotal time:  " << flush;
    MyTime(t.user_time());
    std::cout << std::endl;
  }
}
