#include <iostream>

#include "tests/catch.hpp"

// qvm matrix headers
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include "../../include/ANTL/XGCD/xgcd_plain.hpp"
#include "../../include/ANTL/XGCD/xgcd_binary_l2r.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <vector>
#include <cstdlib>

using namespace NTL;

using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using std::cout;
using std::endl;
using std::vector;
NTL_CLIENT

void set_rand_nbit_value(const int n, int64_t & x){
    int64_t lb = pow(2, n-1);
    int64_t ub = pow(2,n);
    int sign = rand()%2;
    x = (pow(-1, sign))* ((rand()%(ub-lb)) +lb); 
}

/**TEST_CASE("rand nums of bitsize n"){
    #define RAND_MAX pow(2,64);
    srand(time(NULL));
    int64_t x;
    for (int n=5;n<64;n++){
        cout << "n=" << n;
        set_rand_nbit_value(n,x);
        cout << " x=" << x << endl;
        bool positive_case = (x>=pow(2,n-1)) && (x < pow(2,n));
        bool negative_case = (x>-pow(2,n)) && (x<=-pow(2,n-1));
        REQUIRE((positive_case || negative_case));
    }
}
*/

TEST_CASE("XGCD_PLAIN int64_t Benchmarks", "[XGCD][XGCD_PLAIN][Benchmark][int64_t]"){

    #define RAND_MAX pow(2,64);
    srand(time(NULL));
    int64_t g,x,y,a,b;
    int runs = 100;
    int i;
    int numbits;
    BENCHMARK_ADVANCED("XGCD_PLAIN Random 100 run 5-bit inputs")(Catch::Benchmark::Chronometer meter){
        i=0;
        numbits = 5;
        while(i<runs){
            set_rand_nbit_value(numbits, a);
            set_rand_nbit_value(numbits, b);
            meter.measure([] { return XGCD_PLAIN(g,x,y,a,b)});
            i++;
        }
    };
}