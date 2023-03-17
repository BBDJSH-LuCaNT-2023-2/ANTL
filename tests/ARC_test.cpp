#ifndef ARC_TEST
#define ARC_TEST

#define FILE_INDEX 0
#define PROCESS 1
#define NO_PROCESS 2
#define FILE_DONE 3
#define SLAVE_DONE 4

// #define LIST_SIZE_QUADRATIC 33554432
#define LIST_SIZE_QUADRATIC 10000

#include "Quadratic/Applications/Tabulation/iq-data.hpp"
#include "Quadratic/Applications/Tabulation/iq-stats.hpp"
#include <ANTL/Quadratic/QuadraticOrder.hpp>
#include <dirent.h>
#include <fstream>
#include <list>
#include <unistd.h>

#include "Quadratic/Applications/Tabulation/AuxillaryFunctions.hpp"

using namespace NTL;

bool DBG_LENSTRA_TEST = true;


int main (int argc, char ** argv) {

  std::cout << sizeof(long long) << std::endl;
  std::cout << sizeof(NTL::ZZ) << std::endl;

  iq_data combined_data;

  return 0;
}

#endif
