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
#define LIST_SIZE_QUADRATIC 100000000

#include "AuxillaryFunctions.hpp"
#include "timer.hpp"

#include <string>

#include <ANTL/Quadratic/QuadraticOrder.hpp>

NTL_CLIENT
using namespace ANTL;
using namespace std::chrono;

int main(int argc, char **argv) {

  ZZ* Dlist = new ZZ[LIST_SIZE_QUADRATIC];
//   ZZ Dlist[LIST_SIZE_QUADRATIC];
//   std::vector<ZZ> Dlist; Dlist.resize(LIST_SIZE_QUADRATIC);
  ZZ L, H, D, pmax;
  vec_ZZ Cl;
  //   QuadraticOrder<ZZ> QO;
  //   QuadraticOrder<long> QOl;
  //   QuadraticOrder<long long> QOll;
  long long Dll;
  long n, i, Dl;
  //   timer t;
  int vb, tp, table, uncond, alg, prec;

  if (argc != 9) {
    cerr << "Usage:  qcl_interval L H type table uncond vb alg prec"
         << std::endl;
    cerr << "  L - lower bound for |D|" << std::endl;
    cerr << "  H - upper bound for |D|" << std::endl;
    cerr << "  type = 0 (imag), 1 (real)" << std::endl;
    cerr << "  table - 0 (do not use L(s,X) tables), 1 (use them)" << std::endl;
    cerr << "  uncond - 0 (conditional on ERH), 1 (unconditional)" << std::endl;
    cerr << "  vb - verbosity" << std::endl;
    cerr << "  alg - 0 (bsgs), 1 (bs)" << std::endl;
    cerr << "  prec - 0 (automatic), 1 (long), 2 (long long), >2 (ZZ)"
         << std::endl;
    exit(1);
  }

  conv(L, argv[1]);
  conv(H, argv[2]);
  tp = atoi(argv[3]);

  table = atoi(argv[4]);
  uncond = atoi(argv[5]);
  vb = atoi(argv[6]);
  alg = atoi(argv[7]);
  prec = atoi(argv[8]);

  if (vb > 1) {
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
    if (alg == 0) {
      std::cout << "Using BSGS" << std::endl;
    }
    else {
      std::cout << "Using BS" << std::endl;
    }
  }

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
  //   t.start_timer();

  // generate list of fundamental D with L <= |D| <= H
  if (tp == 0) {
    //     get_Dlist_imag(L, H, Dlist, n, 1); // initialize prime list
    //     get_Dlist_imag(L, H, Dlist, n);
  } else {
    get_Dlist_real(L, H, Dlist, n, 1); // initialize prime list
                                       //     get_Dlist_real(L, H, Dlist, n);
    get_Dlist_real(L, H, Dlist, n, 0);
  }
  if (vb > 1) {
    std::cout << "\n# fields = " << n << std::endl;
  }

  auto start_overall_time = high_resolution_clock::now();
  auto finish_overall_time = high_resolution_clock::now();
  duration<long int, std::ratio<1, 1000000>> dur_overall_time(0);

  auto start_trial_time = high_resolution_clock::now();
  auto finish_trial_time = high_resolution_clock::now();

  duration<long int, std::ratio<1, 1000000>> dur_regulator_trial(0);
  duration<long int, std::ratio<1, 1000000>> dur_class_group_trial(0);

  duration<long int, std::ratio<1, 1000000>> total_dur_regulator_trial(0);
  duration<long int, std::ratio<1, 1000000>> total_dur_class_group_trial(0);

//   std::cout << std::fixed;
  for (i = 0; i < n; ++i) {
    D = Dlist[i];
    //       QO.assign(D);
    if (vb > 0) {
      std::cerr << "Computing case " << i << " of " << n << ": D = " << D << std::endl;
    }

    QuadraticOrder<long> QO{to_long(D)};
    std::cout << D << flush;

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
    l_function.init(to_long(D), 2);
    // ===TEMPORARY INITIALIZATION===

    start_trial_time = high_resolution_clock::now();
    pair<double, ZZ> regulator_and_hstar = get_regulator_and_hstar(QO, l_function);
    finish_trial_time = high_resolution_clock::now();

    dur_regulator_trial = duration_cast<microseconds>(finish_trial_time - start_trial_time);

    double regulator = regulator_and_hstar.first;
    ZZ h_star = regulator_and_hstar.second;

    vector<long> class_group;
    if(alg == 0) {
      // Computing and timing a single class group BSGS computation
      start_trial_time = high_resolution_clock::now();
      class_group = get_class_group_BSGS(QO, regulator, h_star);
      finish_trial_time = high_resolution_clock::now();

      dur_class_group_trial = duration_cast<microseconds>(finish_trial_time - start_trial_time);
    }
    else if(alg == 1) {
      // Computing and timing a single class group BS computation
      start_trial_time = high_resolution_clock::now();
      class_group = get_class_group_BS(QO, regulator, h_star);
      finish_trial_time = high_resolution_clock::now();

      dur_class_group_trial = duration_cast<microseconds>(finish_trial_time - start_trial_time);
    }

    //Formatting the class group first
    class_group.erase(std::remove(class_group.begin(), class_group.end(), 1), class_group.end());

    //Computing new durations totals
    total_dur_regulator_trial += dur_regulator_trial;
    total_dur_class_group_trial += dur_class_group_trial;

    //Outputting to stream
    std::cout << " " << long(floor(regulator*1000)) << flush;
    std::cout << " " << class_group << std::endl;
  }

  finish_overall_time = high_resolution_clock::now();
  dur_overall_time = duration_cast<microseconds>(finish_overall_time - start_overall_time);

  // Converting duration to Mins_Secs
  long reg_mins, clg_mins, tot_mins;
  double reg_secs, clg_secs, tot_secs;
  std::string reg_str, clg_str, tot_str;

  reg_mins = total_dur_regulator_trial.count() / 60000000;
  clg_mins = total_dur_class_group_trial.count() / 60000000;
  tot_mins = dur_overall_time.count() / 60000000;

  reg_secs = double(double(total_dur_regulator_trial.count()) / double(1000000)) - double(reg_mins*60);
  clg_secs = double(double(total_dur_class_group_trial.count()) / double(1000000)) - double(clg_mins*60);
  tot_secs = double(double(dur_overall_time.count()) / double(1000000)) - double(tot_mins*60);

  reg_str = to_string(reg_mins);
  clg_str = to_string(clg_mins);
  tot_str = to_string(tot_mins);

  reg_str += "m" + to_string(long(floor(reg_secs))) + ".";
  clg_str += "m" + to_string(long(floor(clg_secs))) + ".";
  tot_str += "m" + to_string(long(floor(tot_secs))) + ".";

  std::string reg_temp, clg_temp, tot_temp;

  reg_temp = to_string(long(floor((reg_secs - floor(reg_secs))*1000))) + "s";
  clg_temp = to_string(long(floor((clg_secs - floor(clg_secs))*1000))) + "s";
  tot_temp = to_string(long(floor((tot_secs - floor(tot_secs))*1000))) + "s";

  reg_temp.insert(reg_temp.begin(), 4 - reg_temp.size(), '0');
  clg_temp.insert(clg_temp.begin(), 4 - clg_temp.size(), '0');
  tot_temp.insert(tot_temp.begin(), 4 - tot_temp.size(), '0');

  reg_str += reg_temp;
  clg_str += clg_temp;
  tot_str += tot_temp;

//   std::cout << "Time spent computing regulator: " << total_dur_regulator_trial.count() / 1000 << std::endl;
//   std::cout << "Time spent computing class grp: " << total_dur_class_group_trial.count() / 1000 << std::endl;
//   std::cout << "Time spent computing overall  : " << dur_overall_time.count() / 1000 << std::endl;
//
//   std::cout << std::endl;

  std::cout << "Time spent computing regulator: " << reg_str << std::endl;
  std::cout << "Time spent computing class grp: " << clg_str << std::endl;
  std::cout << "Time spent computing overall  : " << tot_str << std::endl;


}
