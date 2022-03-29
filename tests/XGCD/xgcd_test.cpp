#include <iostream>

#include "tests/catch.hpp"

// qvm matrix headers
#include <boost/qvm/mat.hpp>
#include <boost/qvm/mat_traits.hpp>
#include <boost/qvm/mat_access.hpp>
#include <boost/qvm/mat_operations.hpp>
#include "../../include/ANTL/XGCD/xgcd_plain.hpp"
#include <boost/math/tools/polynomial.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <vector>

using namespace NTL;

using namespace boost::multiprecision;
using boost::math::tools::polynomial;
using boost::multiprecision::mpf_float;
using std::cout;
using std::endl;
using std::vector;
NTL_CLIENT


/**
 * @brief Abstract Class for testing XGCD
 * @author Zack Baker
 * 
 * 
 * @param g the gcd of a and b
 * @param a first input parameter
 * @param b second input parameter
 * @param u Bezout coefficient for a
 * @param v Bezout coefficient for b
 * 
 * @tparam T The datatype to test 
 * 
 * Abstract class to form the foundation of a class-based testing framework for XGCD algorithms.
 * Generically, all XGCD methods require 5 parameters, intended as follows:
 */
template<typename T>
class XGCDTestInstance {
    protected:
        T g;
        T u;
        T v;
        T a;
        T b; 
    public:
        /**
         * @brief Parent constructor, should be sufficient for all implementing classes
         * 
         * @param a input parameter a
         * @param b input parameter b
         */
        XGCDTestInstance(T a, T b){
            this->a = a;
            this->b = b;
        }
        /**
         * @brief evaluation method for this XGCD class type
         * 
         * Implementing classes must implement this method to evaluate their XGCD algorithm,
         * setting members equal to the evaluation results as appropriate
         * 
         */
        virtual void evaluateXGCD();
        /**
         * @brief Generic test method. Compares the input vector to the stored vector of class members 
         * 
         * @param expected input vector of expected values
         * @return true the expected members all match the actual members
         * @return false the expected vector differed from the member vector
         */
        bool testXGCD(vector<T> expected){
            return expected == createMemberVector(); 
        }
        /**
         * @brief Create a vector of member variables
         * 
         * @return A vector of member variables, in order (g,u,a,v,b)  
         */
        vector<T> createMemberVector(){
            vector<T> members;
            members.push_back(g);
            members.push_back(u);
            members.push_back(a);
            members.push_back(v);
            members.push_back(b);
            return members;
        }

        /**
         * @brief reset the values of the test object, providing new input values a and b
         * 
         * @note Reimplement this method if template assignment to 0 doesn't make sense
         * 
         * @param a first new input 
         * @param b second new input
         */
        void refreshInstance(T a, T b){
            g = 0;
            u = 0;
            v = 0;
            this->a = a;
            this->b = b;
            
        }
         
};

/**
 * @brief Class template for testing XGCDPlain 
 * @author Zack Baker
 * 
 * 
 * @tparam T The datatype to test
 */
template<typename T>
class XGCDPlainTestInstance: public XGCDTestInstance<T> {
    public:
        /**
         * @brief Construct a new XGCDPlainTestInstance object
         * 
         * @param a first input parameter
         * @param b second input parameter
         * 
         * Inherits parent constructor
         */
        XGCDPlainTestInstance(T a, T b) : XGCDTestInstance<T>(a,b){}
        /**
         * @brief Evaluate this XGCDPlainTestInstance by calling XGCD_PLAIN
         * 
         */
        void evaluateXGCD(){
            //hooray for templating
            XGCD_PLAIN(this->g,this->u,this->v,this->a,this->b);
        }
        /**
         * @brief Simplified test method
         * 
         * @param gcd the expected GCD method
         * @return true g is the expected gcd of a and b 
         * @return false g is not the expected gcd of a and b
         */
        bool testXGCD(T gcd){
            evaluateXGCD();
            return gcd == this->g;
        }

};

/*int main(){
long a,b,r,s,g;
a = 3;
b= 5;
XGCD(g, r,s,a,b);
cout << g << " = " << a << "*" << r <<  " + " << b << "*" << s << endl;

XGCD_PLAIN(r,s,a,b);
cout << r << " " << s << " " << a << " " << b << "\n";
return 0;
}*/

TEST_CASE("[XGCD]: int64_t PLAIN"){
    int64_t a,b;

    //basic initial test: (3,5) = 1
    a = 3;
    b = 5;

    //XGCDTestInstance<int>* inst = new XGCDTestInstance<int>(3,4);
    XGCDPlainTestInstance<int64_t>* inst = new XGCDPlainTestInstance<int64_t>(3,5);
    REQUIRE(inst->testXGCD(1));

    inst->refreshInstance(5,3);
    REQUIRE(inst->testXGCD(1));

    //common factor test: (5,20) = 5
    inst->refreshInstance(5,20);
    REQUIRE(inst->testXGCD(5));

    //unit test: (1, 999) = 1
    inst->refreshInstance(1,999);
    REQUIRE(inst->testXGCD(1));

    //sign permutation tests: (+/- 3, +/- 6) = 3
    inst->refreshInstance(-3,6);
    REQUIRE(inst->testXGCD(3));

    inst->refreshInstance(3,-6);
    REQUIRE(inst->testXGCD(3));

    inst->refreshInstance(-3,-6);
    REQUIRE(inst->testXGCD(3));

    //arbitrary large (2^32 < a,b < 2^63) positive test: (3166167471260038366, 2078992898117306689) = 1
    inst->refreshInstance(3166167471260038366, 2078992898117306689);
    REQUIRE(inst->testXGCD(1));

    //arbitrary large (-2^32 > a,b > 2^63) negative test: (-3867470587490682194, -6531477986582055176) = 2 
    inst->refreshInstance(-3867470587490682194, -6531477986582055176);
    REQUIRE(inst->testXGCD(2));

    //0 value tests
    inst->refreshInstance(0,100);
    REQUIRE(inst->testXGCD(100));

    inst->refreshInstance(50,0);
    REQUIRE(inst->testXGCD(50));

    inst->refreshInstance(0,0);
    REQUIRE(inst->testXGCD(0));

}
